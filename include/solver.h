#ifndef SOLVER_H
#define SOLVER_H

#ifdef __cplusplus
extern "C" {
#endif

// Get Ceres Solver major version number
int test_ceres(void);

// Run Ceres hello world example (minimizes 0.5 * (10 - x)^2)
double run_hello_world(void);

#ifdef __cplusplus
}
#endif

#endif // SOLVER_H