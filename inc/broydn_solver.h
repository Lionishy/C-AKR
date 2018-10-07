#ifndef IKI_BroydnSolver_H
#define IKI_BroydnSolver_H

#include <solver.h>

IterativeSolverStatus broydn_solver(
    func_2d F, void const * Fcntx, double const * initial_guess
    , double dx, double dy
    , double eps, unsigned int max_iteration
    , double * res_vec, unsigned int * res_iteration
);

#endif //IKI_BroydnSolver_H