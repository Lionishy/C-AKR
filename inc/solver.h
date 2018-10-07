#ifndef IKI_Solver_H
#define IKI_Solver_H

#include <complex.h>

typedef void (*func)(void const *context, double in, double *out);
typedef void (*func_2d)(void const *context, double const * in_vec, double * out_vec);
typedef void (*func_cpx)(void const *context, double complex in, double complex *out);

typedef enum _IterativeSolverStatus {
      SLV_OK = 0
    , SLV_MAX_ITER
    , SLV_NO_PROGRESS
    , SLV_SINGULARITY
} IterativeSolverStatus;

typedef IterativeSolverStatus (*solver) (
    void const *context_ptr
    , func F, void const *F_context, double initial_guess
    , double *res, unsigned *res_iteration
);

typedef IterativeSolverStatus (*solver_2d) (
      void const *env_ptr
    , func_2d F, void const * F_context, double const * initial_guess
    , double * res_vec, unsigned int * res_iteration
);

#endif