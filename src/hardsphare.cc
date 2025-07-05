#include <ceres.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <unsupported/Eigen/FFT>
#include <vector>

class HardSphereDFT : public ceres::CostFunction {
   public:
    struct Config {
        double sigma;               // 硬球直径
        double temperature;         // 温度
        double chemical_potential;  // 化学ポテンシャル
        double dr;                  // 格子間隔
        int num_points;             // 格子点数
        double r_max;               // 最大距離
    };

    explicit HardSphereDFT(const Config& config)
        : config_(config), fft_initialized_(false) {
        // Ceresの設定
        set_num_residuals(config.num_points);
        mutable_parameter_block_sizes()->push_back(config.num_points);

        // 格子の初期化
        initialize_grid();

        // 重み関数の計算
        compute_weight_functions();

        // FFTの初期化
        initialize_fft();

        // Jacobian行列の事前計算
        precompute_jacobian_matrices();
    }

    ~HardSphereDFT() { cleanup_fft(); }

    bool Evaluate(double const* const* parameters,
                  double* residuals,
                  double** jacobians) const override {
        const double* rho = parameters[0];
        const int n = config_.num_points;

        // 密度の物理的制約チェック
        for (int i = 0; i < n; ++i) {
            if (rho[i] < 0.0) {
                // 負密度のペナルティ
                for (int j = 0; j < n; ++j) {
                    residuals[j] = (j == i) ? -rho[i] * 1e6 : 0.0;
                }
                if (jacobians && jacobians[0]) {
                    std::fill(jacobians[0], jacobians[0] + n * n, 0.0);
                    jacobians[0][i * n + i] = -1e6;
                }
                return true;
            }
        }

        // 関数値計算：FFTで重み密度を計算
        compute_residuals_with_fft(rho, residuals);

        // Jacobian計算：事前計算された行列を使用
        if (jacobians && jacobians[0]) {
            compute_jacobian(rho, jacobians[0]);
        }

        return true;
    }

   private:
    Config config_;
    std::vector<double> r_grid_;

    // 重み関数（Rosenfeld functional用）
    std::vector<double> w0_, w1_, w2_, w3_;

    // FFT関連
    mutable bool fft_initialized_;
    mutable Eigen::FFT<double> fft_;
    mutable std::vector<std::complex<double>> fft_w0_, fft_w1_, fft_w2_,
        fft_w3_;
    mutable std::vector<std::complex<double>> fft_rho_, fft_result_;
    mutable int padded_size_;

    // 事前計算されたJacobian成分
    std::vector<std::vector<double>> jacobian_w0_, jacobian_w1_, jacobian_w2_,
        jacobian_w3_;

    void initialize_grid() {
        const int n = config_.num_points;
        r_grid_.resize(n);

        for (int i = 0; i < n; ++i) {
            r_grid_[i] = (i + 0.5) * config_.dr;  // セル中心での評価
        }
    }

    void compute_weight_functions() {
        const int n = config_.num_points;
        const double sigma = config_.sigma;
        const double R = sigma / 2.0;  // 硬球半径

        w0_.resize(n);
        w1_.resize(n);
        w2_.resize(n);
        w3_.resize(n);

        for (int i = 0; i < n; ++i) {
            const double r = r_grid_[i];

            if (r <= R) {
                // Rosenfeld重み関数
                w3_[i] = 1.0;
                w2_[i] = 1.0;
                w1_[i] = 1.0 / (4.0 * M_PI * R);
                w0_[i] = 1.0 / (4.0 * M_PI * R * R);
            } else {
                w0_[i] = w1_[i] = w2_[i] = w3_[i] = 0.0;
            }
        }
    }

    void initialize_fft() const {
        if (fft_initialized_) return;

        const int n = config_.num_points;
        padded_size_ = next_power_of_2(2 * n);

        // メモリ確保
        fft_w0_.resize(padded_size_);
        fft_w1_.resize(padded_size_);
        fft_w2_.resize(padded_size_);
        fft_w3_.resize(padded_size_);
        fft_rho_.resize(padded_size_);
        fft_result_.resize(padded_size_);

        // Eigen FFTはプラン不要

        // 重み関数のFFT（一度だけ計算）
        fft_weight_function(w0_, fft_w0_);
        fft_weight_function(w1_, fft_w1_);
        fft_weight_function(w2_, fft_w2_);
        fft_weight_function(w3_, fft_w3_);

        fft_initialized_ = true;
    }

    void fft_weight_function(
        const std::vector<double>& weight,
        std::vector<std::complex<double>>& fft_weight) const {
        const int n = config_.num_points;

        // ゼロパディング
        std::fill(fft_weight.begin(),
                  fft_weight.end(),
                  std::complex<double>(0.0, 0.0));

        // 重み関数を体積要素込みで設定
        for (int i = 0; i < n; ++i) {
            double volume_element =
                4.0 * M_PI * r_grid_[i] * r_grid_[i] * config_.dr;
            double simpson_weight = compute_simpson_weight(i, n) / 3.0;
            fft_weight[i] = std::complex<double>(
                weight[i] * volume_element * simpson_weight, 0.0);
        }

        fft_.fwd(fft_weight, fft_weight);
    }

    void compute_residuals_with_fft(const double* rho,
                                    double* residuals) const {
        const int n = config_.num_points;

        // 密度のFFT
        compute_density_fft(rho);

        // 各重み密度の畳み込み計算
        std::vector<double> n0(n), n1(n), n2(n), n3(n);

        compute_weighted_density(fft_w0_, n0);
        compute_weighted_density(fft_w1_, n1);
        compute_weighted_density(fft_w2_, n2);
        compute_weighted_density(fft_w3_, n3);

        // Rosenfeld汎関数による化学ポテンシャル計算
        for (int i = 0; i < n; ++i) {
            // 理想気体部分
            double mu_ideal = std::log(rho[i]);

            // 過剰部分（Rosenfeld functional）
            double xi = 1.0 - n3[i];
            double xi2 = xi * xi;
            double xi3 = xi2 * xi;

            if (xi <= 0.0) {
                // 密度が高すぎる場合のペナルティ
                residuals[i] = 1e6;
                continue;
            }

            double mu_excess = -std::log(xi) + n2[i] / xi +
                               n1[i] * n2[i] / xi2 +
                               n0[i] * n2[i] * n2[i] / xi3;

            // 変分条件：μ_total = μ_chemical
            residuals[i] = mu_ideal + mu_excess - config_.chemical_potential;
        }
    }

    void compute_density_fft(const double* rho) const {
        const int n = config_.num_points;

        // ゼロパディング
        std::fill(
            fft_rho_.begin(), fft_rho_.end(), std::complex<double>(0.0, 0.0));

        // 密度を設定
        for (int i = 0; i < n; ++i) {
            fft_rho_[i] = std::complex<double>(rho[i], 0.0);
        }

        fft_.fwd(fft_rho_, fft_rho_);
    }

    void compute_weighted_density(
        const std::vector<std::complex<double>>& fft_weight,
        std::vector<double>& weighted_density) const {
        const int n = config_.num_points;

        // 周波数域での乗算
        for (int i = 0; i < padded_size_; ++i) {
            fft_result_[i] = fft_rho_[i] * fft_weight[i];
        }

        // 逆FFT
        fft_.inv(fft_result_, fft_result_);

        // 結果の抽出と正規化
        weighted_density.resize(n);
        for (int i = 0; i < n; ++i) {
            weighted_density[i] = fft_result_[i].real();
        }
    }

    void precompute_jacobian_matrices() {
        const int n = config_.num_points;

        jacobian_w0_.resize(n, std::vector<double>(n));
        jacobian_w1_.resize(n, std::vector<double>(n));
        jacobian_w2_.resize(n, std::vector<double>(n));
        jacobian_w3_.resize(n, std::vector<double>(n));

        // ∂n_α[i]/∂ρ[k] = w_α(|r_i - r_k|) * volume_element[k] *
        // integration_weight[k]
        for (int i = 0; i < n; ++i) {
            for (int k = 0; k < n; ++k) {
                double r_ik = std::abs(r_grid_[i] - r_grid_[k]);

                double volume_element =
                    4.0 * M_PI * r_grid_[k] * r_grid_[k] * config_.dr;
                double simpson_weight = compute_simpson_weight(k, n) / 3.0;
                double total_weight = volume_element * simpson_weight;

                jacobian_w0_[i][k] =
                    interpolate_weight_function(w0_, r_ik) * total_weight;
                jacobian_w1_[i][k] =
                    interpolate_weight_function(w1_, r_ik) * total_weight;
                jacobian_w2_[i][k] =
                    interpolate_weight_function(w2_, r_ik) * total_weight;
                jacobian_w3_[i][k] =
                    interpolate_weight_function(w3_, r_ik) * total_weight;
            }
        }
    }

    void compute_jacobian(const double* rho, double* jacobian) const {
        const int n = config_.num_points;

        // 現在の重み密度を計算（FFT使用）
        std::vector<double> n0(n), n1(n), n2(n), n3(n);
        compute_density_fft(rho);
        compute_weighted_density(fft_w0_, n0);
        compute_weighted_density(fft_w1_, n1);
        compute_weighted_density(fft_w2_, n2);
        compute_weighted_density(fft_w3_, n3);

        // Jacobian計算
        for (int i = 0; i < n; ++i) {
            for (int k = 0; k < n; ++k) {
                double jac = 0.0;

                if (i == k) {
                    // 理想気体部分の対角項
                    jac += 1.0 / rho[i];
                }

                // 過剰部分の微分
                double xi = 1.0 - n3[i];
                if (xi > 0.0) {
                    double xi2 = xi * xi;
                    double xi3 = xi2 * xi;
                    double xi4 = xi3 * xi;

                    // ∂μ_excess[i]/∂n_α[i] の計算
                    double dmu_dn0 = n2[i] * n2[i] / xi3;
                    double dmu_dn1 = n2[i] / xi2;
                    double dmu_dn2 =
                        1.0 / xi + n1[i] / xi2 + 2.0 * n0[i] * n2[i] / xi3;
                    double dmu_dn3 = 1.0 / xi + n2[i] / xi2 +
                                     2.0 * n1[i] * n2[i] / xi3 +
                                     3.0 * n0[i] * n2[i] * n2[i] / xi4;

                    // チェインルール適用
                    jac += dmu_dn0 * jacobian_w0_[i][k] +
                           dmu_dn1 * jacobian_w1_[i][k] +
                           dmu_dn2 * jacobian_w2_[i][k] +
                           dmu_dn3 * jacobian_w3_[i][k];
                }

                jacobian[i * n + k] = jac;
            }
        }
    }

    double interpolate_weight_function(const std::vector<double>& weight,
                                       double r) const {
        if (r >= config_.sigma / 2.0) return 0.0;

        // 最も近い格子点での値を返す（簡単な実装）
        int index = static_cast<int>(r / config_.dr);
        if (index >= weight.size()) return 0.0;
        return weight[index];
    }

    double compute_simpson_weight(int j, int n) const {
        if (j == 0 || j == n - 1) {
            return 1.0;
        } else if (j % 2 == 1) {
            return 4.0;
        } else {
            return 2.0;
        }
    }

    int next_power_of_2(int n) const {
        int power = 1;
        while (power < n) power *= 2;
        return power;
    }

    void cleanup_fft() const {
        if (!fft_initialized_) return;
        // Eigen FFTは自動的にメモリ管理されるのでクリーンアップ不要
    }
};

// 使用例
class HardSphereDFTSolver {
   public:
    static std::vector<double> solve(
        const HardSphereDFT::Config& config,
        const std::vector<double>& initial_density) {
        std::vector<double> rho = initial_density;

        // Ceresの問題設定
        ceres::Problem problem;

        auto cost_function = new HardSphereDFT(config);
        problem.AddResidualBlock(cost_function, nullptr, rho.data());

        // 密度の正値制約
        for (int i = 0; i < config.num_points; ++i) {
            problem.SetParameterLowerBound(rho.data(), i, 1e-12);
        }

        // ソルバー設定
        ceres::Solver::Options options;
        options.linear_solver_type = ceres::DENSE_QR;
        options.minimizer_progress_to_stdout = true;
        options.max_num_iterations = 100;
        options.function_tolerance = 1e-10;
        options.gradient_tolerance = 1e-12;

        ceres::Solver::Summary summary;
        ceres::Solve(options, &problem, &summary);

        std::cout << summary.BriefReport() << std::endl;

        return rho;
    }
};

// 簡単な使用例
/*
int main() {
    // 硬球系の設定
    HardSphereDFT::Config config;
    config.sigma = 1.0;              // 硬球直径
    config.temperature = 1.0;        // 温度
    config.chemical_potential = 2.0; // 化学ポテンシャル
    config.dr = 0.05;               // 格子間隔
    config.num_points = 512;        // 格子点数
    config.r_max = config.dr * config.num_points;

    // 初期密度（均一密度で開始）
    std::vector<double> initial_density(config.num_points, 0.5);

    // DFT計算実行
    auto solution = HardSphereDFTSolver::solve(config, initial_density);

    // 結果の出力
    std::cout << "Final density profile:" << std::endl;
    for (int i = 0; i < std::min(10, (int)solution.size()); ++i) {
        double r = (i + 0.5) * config.dr;
        std::cout << "r=" << r << ", rho=" << solution[i] << std::endl;
    }

    return 0;
}
*/
