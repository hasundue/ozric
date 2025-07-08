#include "ceres.h"

#include <ceres_wrapped.h>

#include <cmath>
#include <vector>

class ThreadLocal {
   public:
    void* data;

    ThreadLocal(thread_local_init_t init, thread_local_deinit_t deinit)
        : deinit_(deinit) {
        data = init();
    }
    ~ThreadLocal() { deinit_(data); }

   private:
    thread_local_deinit_t deinit_;
};

class OzricCostFunction : public ceres::CostFunction {
   public:
    OzricCostFunction(cost_function_t cost_function,
                      thread_local_init_t thread_local_init,
                      thread_local_deinit_t thread_local_deinit,
                      int num_points)
        : cost_function_(cost_function),
          thread_local_init_(thread_local_init),
          thread_local_deinit_(thread_local_deinit) {
        set_num_residuals(num_points);
        mutable_parameter_block_sizes()->push_back(num_points);
    }
    ~OzricCostFunction() {}

    bool Evaluate(double const* const* parameters,
                  double* residuals,
                  double** jacobians) const override {
        thread_local ThreadLocal tl(thread_local_init_, thread_local_deinit_);
        return (*cost_function_)(
            tl.data, const_cast<double**>(parameters), residuals, jacobians);
    }

   private:
    cost_function_t cost_function_;
    thread_local_init_t thread_local_init_;
    thread_local_deinit_t thread_local_deinit_;
};

extern "C" void* ceres_create_cost_function(
    cost_function_t cost_function,
    thread_local_init_t thread_local_init,
    thread_local_deinit_t thread_local_deinit,
    int num_points) {
    return new OzricCostFunction(
        cost_function, thread_local_init, thread_local_deinit, num_points);
}

extern "C" void* ceres_create_problem() { return new ceres::Problem(); }

extern "C" void ceres_problem_add_residual_block(void* problem,
                                                 void* cost_function,
                                                 double** parameters,
                                                 int num_parameter_blocks) {
    ceres::Problem* p = static_cast<ceres::Problem*>(problem);
    ceres::CostFunction* cf = static_cast<ceres::CostFunction*>(cost_function);

    std::vector<double*> parameter_blocks;
    for (int i = 0; i < num_parameter_blocks; i++) {
        parameter_blocks.push_back(parameters[i]);
    }

    p->AddResidualBlock(cf, nullptr, parameter_blocks);
}

// Solve the optimization problem
extern "C" int ceres_problem_solve(void* problem) {
    ceres::Problem* p = static_cast<ceres::Problem*>(problem);

    // Set up solver options
    ceres::Solver::Options options;
    options.linear_solver_type = ceres::DENSE_QR;
    options.minimizer_progress_to_stdout = false;  // Disable verbose output

    // Solve the problem
    ceres::Solver::Summary summary;
    ceres::Solve(options, p, &summary);

    // Return 1 for success, 0 for failure
    return summary.termination_type == ceres::CONVERGENCE ? 1 : 0;
}

// Hello World cost functor for auto-differentiation
struct HelloWorldCostFunctor {
    template <typename T>
    bool operator()(const T* const x, T* residual) const {
        residual[0] = T(10.0) - x[0];
        return true;
    }
};

// Export instance of HelloWorldCostFunction using auto-differentiation
extern "C" void* ceres_create_hello_world_cost_function() {
    return new ceres::AutoDiffCostFunction<HelloWorldCostFunctor, 1, 1>(
        new HelloWorldCostFunctor);
}
