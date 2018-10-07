#ifndef IKI_SecantSolver_H
#define IKI_SecantSolver_H

#include <solver.h>

IterativeSolverStatus secant_solver(
      func F, void const *Fcntx, double initial_guess
    , double dx, double eps, unsigned max_iter
    , double *res, unsigned *iterations
);

#endif //IKI_SecantSolver_H