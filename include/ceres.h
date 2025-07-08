#pragma once

#ifdef __cplusplus
extern "C" {
#endif

// Function pointer types for cost function callbacks
typedef void* (*thread_local_init_t)();
typedef void (*thread_local_deinit_t)(void*);

typedef int (*cost_function_t)(void* data,
                               double** parameters,
                               double* residuals,
                               double** jacobians);

typedef cost_function_t (*create_cost_function_t)(void* point);

// Create OzricCostFunction and return as opaque pointer
void* ceres_create_cost_function(cost_function_t cost_function,
                                 thread_local_init_t thread_local_init,
                                 thread_local_deinit_t thread_local_deinit,
                                 int num_points);

// Create ceres::Problem and return as opaque pointer
void* ceres_create_problem();

// Add residual block to ceres::Problem
void ceres_problem_add_residual_block(void* problem,
                                      void* cost_function,
                                      double** parameter_blocks,
                                      int num_parameter_blocks);

// Solve the optimization problem (returns 1 for success, 0 for failure)
int ceres_problem_solve(void* problem);

// Create HelloWorldCostFunction instance and return as opaque pointer
void* ceres_create_hello_world_cost_function();

#ifdef __cplusplus
}
#endif
