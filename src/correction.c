#include <correction.h>
#include <vector.h>
#include <physical_environment.h>
#include <broydn_solver.h>
#include <secant_solver.h>

#include <stddef.h>
#include <stdbool.h>
#include <math.h>

bool trivial_corrector(void const *cntx, VectorSp *R, VectorSp *K, double *w_vec) {
    return true;
}

typedef struct _omega_closureH_context {
    VectorSp R;
    VectorH Kh;
    relationH_t relation;
    void const *relation_cntx;
} omega_closureH_context_t;

static inline
void omega_closureH(void const *ptr, double const *in_vec, double *out_vec) {
    omega_closureH_context_t const *cntx = ptr;
    cntx->relation(cntx->relation_cntx,cntx->R,cntx->Kh,in_vec,out_vec);
}

bool omega_correctorH(void const *ptr, VectorSp *R, VectorH *K, double *w_vec) {
    omega_correctorH_context_t const *cntx = ptr;
    omega_closureH_context_t closure_cntx = {
          *R, *K
        , cntx->relation, cntx->relation_cntx
    };

    double initial_guess[2] = {w_vec[0],w_vec[1]}, res_vec[2]; unsigned iterations;
    IterativeSolverStatus res_status =  broydn_solver(
            omega_closureH, &closure_cntx, initial_guess
            ,cntx->dx, cntx->dy, cntx->eps
            ,cntx->max_iter
            ,res_vec, &iterations
        );

    if (SLV_OK == res_status) {
        w_vec[0] = res_vec[0]; w_vec[1] = res_vec[1]; 
        return true;
    }

    return false;
}

typedef struct _omega_closureSp_context {
    VectorSp R;
    VectorSp K;
    relationSp_t relation;
    void const *relation_cntx;
} omega_closureSp_context_t;

static inline
void omega_closureSp(void const *ptr, double const *in_vec, double *out_vec) {
    omega_closureSp_context_t const *cntx = ptr;
    cntx->relation(cntx->relation_cntx,cntx->R,cntx->K,in_vec,out_vec);
}

bool omega_correctorSp(void const *ptr, VectorSp *R, VectorSp *K, double *w_vec) {
    omega_correctorSp_context_t const *cntx = ptr;
    omega_closureSp_context_t closure_cntx = {
          *R, *K
        , cntx->relation, cntx->relation_cntx
    };

    double initial_guess[2] = {w_vec[0],w_vec[1]}, res_vec[2]; unsigned iterations;
    IterativeSolverStatus res_status =  broydn_solver(
            omega_closureSp, &closure_cntx, initial_guess
            ,cntx->dx, cntx->dy, cntx->eps
            ,cntx->max_iter
            ,res_vec, &iterations
        );

    if (SLV_OK == res_status) {
        w_vec[0] = res_vec[0]; w_vec[1] = res_vec[1]; 
        return true;
    }

    return false;
}


typedef struct _gamma_kpr_closureH_context {
    VectorSp R;
    VectorH Kh;
    double w_vec[2];
    relationH_t relation;
    void const *relation_cntx;
} gamma_kpr_closureH_context_t;

//in_vec[0]*K0.pr = Kpr  in_vec[1]*w_vec[1] = gamma
static inline
void gamma_kpr_closureH(void const * ptr, double const * in_vec, double * out_vec) {
    gamma_kpr_closureH_context_t const *cntx = ptr;

    double w_vec[2]; w_vec[0] = cntx->w_vec[0]; w_vec[1] = cntx->w_vec[1]*in_vec[1];
    cntx->relation(cntx->relation_cntx,cntx->R,(VectorH){cntx->Kh.pl,cntx->Kh.pr*in_vec[0]},w_vec,out_vec);
}

bool gamma_kpr_correctorH(void const *ptr, VectorSp *R, VectorH *K, double *w_vec) {
    gamma_kpr_correctorH_context_t const *cntx = ptr;
    gamma_kpr_closureH_context_t closure_cntx = {
          *R, *K, {w_vec[0],w_vec[1]}
        , cntx->relation, cntx->relation_cntx
    };

    double initial_guess[2] = {1.,1.}, res_vec[2]; unsigned iterations;
    IterativeSolverStatus res_status =  broydn_solver(
             gamma_kpr_closureH, &closure_cntx, initial_guess
            ,cntx->dx, cntx->dy, cntx->eps
            ,cntx->max_iter
            ,res_vec, &iterations
        );

    if (SLV_OK == res_status) {
        res_vec[0] = res_vec[0] > 0. ? res_vec[0] : -res_vec[0];
        res_vec[1] = res_vec[1] > 0. ? res_vec[1] : -res_vec[1];
        K->pr    *= res_vec[0];
        w_vec[1] *= res_vec[1];

        return true;
    }

    return false;
}

bool gamma_kpr_correctorSp(void const *ptr, VectorSp *R, VectorSp *K, double *w_vec) {
    gamma_kpr_correctorSp_context_t const *cntx = ptr;
    
    PhysicalEnvironment phys_env = cntx->phys_env(cntx->phys_env_cntx,*R);
    VectorSp H = phys_env.H;
    VectorH Kh = projection_of_on(*K,H);
    VectorSp Kpl, Kpr;
    projection_of_on_insp(*K,H,&Kpl,&Kpr);

    gamma_kpr_closureH_context_t closure_cntx = {
          *R, Kh, {w_vec[0],w_vec[1]}
        , cntx->relation, cntx->relation_cntx
    };

    double initial_guess[2] = {1.,1.}, res_vec[2]; unsigned iterations;
    IterativeSolverStatus res_status =  broydn_solver(
             gamma_kpr_closureH, &closure_cntx, initial_guess
            ,cntx->dx, cntx->dy, cntx->eps
            ,cntx->max_iter
            ,res_vec, &iterations
        );

    if (SLV_OK == res_status) {
        w_vec[1] *= res_vec[1];
        *K = (VectorSp){Kpl.r + Kpr.r*res_vec[0], Kpl.th + Kpr.th*res_vec[0], Kpl.phi + Kpr.phi*res_vec[0]};
        return true;
    }

    return false;
}


typedef struct _gamma_kpl_closureH_context {
    VectorSp R;
    VectorH Kh;
    double w_vec[2];
    relationH_t relation;
    void const *relation_cntx;
} gamma_kpl_closureH_context_t;

//in_vec[0]*K0.pr = Kpr  in_vec[1]*w_vec[1] = gamma
static inline
void gamma_kpl_closureH(void const * ptr, double const * in_vec, double * out_vec) {
    gamma_kpl_closureH_context_t const *cntx = ptr;
    
    double w_vec[2]; w_vec[0] = cntx->w_vec[0]; w_vec[1] = cntx->w_vec[1]*in_vec[1];
    cntx->relation(cntx->relation_cntx,cntx->R,(VectorH){cntx->Kh.pl*in_vec[0],cntx->Kh.pr},w_vec,out_vec);
}

bool gamma_kpl_correctorH(void const *ptr, VectorSp *R, VectorH *K, double *w_vec) {
    gamma_kpl_correctorH_context_t const *cntx = ptr;
    
    gamma_kpl_closureH_context_t closure_cntx = {
          *R, *K, {w_vec[0],w_vec[1]}
        , cntx->relation, cntx->relation_cntx
    };

    double initial_guess[2] = {1.,1.}, res_vec[2]; unsigned iterations;
    IterativeSolverStatus res_status =  broydn_solver(
             gamma_kpl_closureH, &closure_cntx, initial_guess
            ,cntx->dx, cntx->dy, cntx->eps
            ,cntx->max_iter
            ,res_vec, &iterations
        );

    if (SLV_OK == res_status) {
        res_vec[0] = res_vec[0] > 0. ? res_vec[0] : -res_vec[0];
        res_vec[1] = res_vec[1] > 0. ? res_vec[1] : -res_vec[1];
        K->pl    *= res_vec[0];
        w_vec[1] *= res_vec[1];

        return true;
    }

    return false;
}

bool gamma_kpl_correctorSp(void const *ptr, VectorSp *R, VectorSp *K, double *w_vec) {
    gamma_kpl_correctorSp_context_t const *cntx = ptr;
    
    PhysicalEnvironment phys_env = cntx->phys_env(cntx->phys_env_cntx,*R);
    VectorSp H = phys_env.H;
    VectorH Kh = projection_of_on(*K,H);
    VectorSp Kpl, Kpr;
    projection_of_on_insp(*K,H,&Kpl,&Kpr);

    gamma_kpl_closureH_context_t closure_cntx = {
          *R, Kh, {w_vec[0],w_vec[1]}
        , cntx->relation, cntx->relation_cntx
    };

    double initial_guess[2] = {1.,1.}, res_vec[2]; unsigned iterations;
    IterativeSolverStatus res_status =  broydn_solver(
             gamma_kpl_closureH, &closure_cntx, initial_guess
            ,cntx->dx, cntx->dy, cntx->eps
            ,cntx->max_iter
            ,res_vec, &iterations
        );

    if (SLV_OK == res_status) {
        w_vec[1] *= res_vec[1];
        *K = (VectorSp){Kpl.r*res_vec[0] + Kpr.r, Kpl.th*res_vec[0] + Kpr.th, Kpl.phi*res_vec[0] + Kpr.phi};
        return true;
    }

    return false;
}


typedef struct _kpl_closureH_context {
    VectorSp R;
    VectorH Kh;
    double w_vec[2];
    relationH_t relation;
    void const *relation_cntx;
} kpl_closureH_context_t;

static inline
void kpl_closureH(void const *ptr, double in, double *out) {
    kpl_closureH_context_t const *cntx = ptr;
    
    double out_vec[2];
    cntx->relation(cntx->relation_cntx,cntx->R,(VectorH){cntx->Kh.pl*in,cntx->Kh.pr},cntx->w_vec,out_vec);
    *out = out_vec[0];
}

bool kpl_correctorH(void const *ptr, VectorSp *R, VectorH *K, double *w_vec) {
    kpl_correctorH_context_t const *cntx = ptr;

    kpl_closureH_context_t closure_cntx = {
         *R, *K
        , {w_vec[0],w_vec[1]}
        , cntx->relation, cntx->relation_cntx
    };

    double initial_guess = 1.; double res; unsigned iterations;
    IterativeSolverStatus res_status =  secant_solver(
          kpl_closureH, &closure_cntx, initial_guess
        , cntx->dx, cntx->eps, cntx->max_iter
        , &res, &iterations
    );

    if (SLV_OK == res_status) {
        K->pl *= res;
        return true;
    }
    
    return false;
}

bool kpl_correctorSp(void const *ptr, VectorSp *R, VectorSp *K, double *w_vec) {
    kpl_correctorSp_context_t const *cntx = ptr;

    PhysicalEnvironment phys_env = cntx->phys_env(cntx->phys_env_cntx,*R);
    VectorSp H = phys_env.H;
    VectorH Kh = projection_of_on(*K,H);
    VectorSp Kpl, Kpr;
    projection_of_on_insp(*K,H,&Kpl,&Kpr);

    kpl_closureH_context_t closure_cntx = {
          *R, Kh, {w_vec[0],w_vec[1]}
        , cntx->relation, cntx->relation_cntx
    };

    double initial_guess = 1.; double res; unsigned iterations;
    IterativeSolverStatus res_status =  secant_solver(
          kpl_closureH, &closure_cntx, initial_guess
        , cntx->dx, cntx->eps, cntx->max_iter
        , &res, &iterations
    );

    if (SLV_OK == res_status) {
        *K = (VectorSp){Kpl.r*res + Kpr.r, Kpl.th*res + Kpr.th, Kpl.phi*res + Kpr.phi};
        return true;
    }

    return false;
}
