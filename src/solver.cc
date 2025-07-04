// Test including ceres headers (now uses our macro-fixed version)
#include <ceres.h>
#include <cmath>

// Simple function that uses ceres header
extern "C" int test_ceres() {
    return CERES_VERSION_MAJOR;
}

// Cost function for hello world example
struct CostFunctor {
    template <typename T>
    bool operator()(const T* const x, T* residual) const {
        residual[0] = T(10.0) - x[0];
        return true;
    }
};

// Real Ceres hello world example - minimize 0.5 * (10 - x)^2
extern "C" double solve() {
    // The variable to solve for, with an initial value.
    double initial_x = 5.0;
    double x = initial_x;

    // Build the problem.
    ceres::Problem problem;

    // Set up the only cost function (also known as residual).
    ceres::CostFunction* cost_function =
        new ceres::AutoDiffCostFunction<CostFunctor, 1, 1>(new CostFunctor);
    problem.AddResidualBlock(cost_function, nullptr, &x);

    // Run the solver!
    ceres::Solver::Options options;
    options.linear_solver_type = ceres::DENSE_QR;
    options.minimizer_progress_to_stdout = false;
    ceres::Solver::Summary summary;
    ceres::Solve(options, &problem, &summary);

    return x;
}
