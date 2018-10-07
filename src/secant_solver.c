#include <secant_solver.h>

#include <math.h>
#include <float.h>

double calculate_derivative(func F, void const *Fcntx, double x0, double dx) {
    double Fp, Fm;
    F(Fcntx,x0+dx,&Fp); F(Fcntx,x0-dx,&Fm);
    return 0.5*(Fp-Fm)/dx;
}

IterativeSolverStatus secant_solver(
      func F, void const *Fcntx, double initial_guess
    , double dx, double eps, unsigned max_iter
    , double *res, unsigned *res_iterations
) {
    IterativeSolverStatus res_status = SLV_OK;
    unsigned int iterations = 0;
    double curr, Fcurr, next, Fnext, G_over;
    curr = initial_guess; F(Fcntx,curr,&Fcurr);
    G_over = 1.0/calculate_derivative(F,Fcntx,curr,dx);

    while (fabs(Fcurr) > eps) {
        if (iterations > max_iter) {
            res_status = SLV_MAX_ITER;
            break;
        }

        if ( isinf(G_over) || isnan(G_over) ) {
            res_status = SLV_SINGULARITY;
            break;
        }

        next = curr - Fcurr*G_over;
        F(Fcntx,next,&Fnext);

        G_over = (next-curr)/(Fnext-Fcurr);
        curr = next;
        Fcurr = Fnext;
        
        ++iterations;
    }

    *res = curr;
    *res_iterations = iterations;
    return res_status;
}