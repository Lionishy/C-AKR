#include <dispersion_relation.h>
#include <vector.h>

#include <complex.h>
#include <math.h>

struct epsilon {
    double complex eps1, eps2;
};

static inline 
struct epsilon source_epsilon(PhysicalEnvironment env, VectorH K, double complex w) {
    double complex w_sqr = w*w;
    double omega_pss = env.omega_ps*env.omega_ps;

    double sf = 1. - 0.5*(env.V.pl*env.V.pl + env.V.pr*env.V.pr);
    double kv_pl = K.pl*env.V.pl;

    double complex
        f1 = w-env.omega_cc*sf-kv_pl,
        f2 = w+env.omega_cc*sf-kv_pl;

    double complex
        F1 = (w-kv_pl)/f1 + 0.5*env.V.pr*env.V.pr*(K.pl*K.pl-w*env.omega_cc)/f1/f1,
        F2 = (w-kv_pl)/f2 + 0.5*env.V.pr*env.V.pr*(K.pl*K.pl-w*env.omega_cc)/f2/f2;

    double complex
        eps1 = -0.5*omega_pss/w_sqr*(F1+F2),
        eps2 =  0.5*omega_pss/w_sqr*(F1-F2);

    return (struct epsilon){eps1,eps2};
}

static inline 
struct epsilon cold_epsilon(PhysicalEnvironment env, VectorH K, double complex w) {
    double complex w_sqr = w*w;
    double omega_pcs = env.omega_pc*env.omega_pc;
    double omega_ccs = env.omega_cc*env.omega_cc;

    double complex
        eps1 = -(omega_pcs)/(w_sqr-omega_ccs),
        eps2 = env.omega_cc/w*omega_pcs/(w_sqr-omega_ccs);

   return (struct epsilon){eps1,eps2}; 
}


static inline
double complex cold_dispersion_relation(PhysicalEnvironment env, VectorH K, double complex w) {
    double complex n_pl = K.pl/w, n_pr = K.pr/w;
    struct epsilon eps_c = cold_epsilon(env,K,w);
    double complex eps1 = 1. + eps_c.eps1, eps2 = eps_c.eps2;
    return (((eps1-n_pl*n_pl)/eps1 - n_pr*n_pr)*((eps1*eps1-eps2*eps2)/eps1 - n_pl*n_pl - n_pr*n_pr) - n_pl*n_pl*eps2/eps1*eps2/eps1);
}

static inline
double complex source_dispersion_relation(PhysicalEnvironment env, VectorH K, double complex w) {
    double complex n_pl = K.pl/w, n_pr = K.pr/w;
    struct epsilon eps_s = source_epsilon(env,K,w);
    double complex eps1 = 1. + eps_s.eps1, eps2 = eps_s.eps2;
    return ((eps1-n_pl*n_pl)/eps1 - n_pr*n_pr)*((eps1*eps1-eps2*eps2)/eps1 - n_pl*n_pl - n_pr*n_pr) - n_pl*n_pl*eps2/eps1*eps2/eps1;
}

static inline
double complex warm_dispersion_relation(PhysicalEnvironment env, VectorH K, double complex w) {
    double complex n_pl = K.pl/w, n_pr = K.pr/w; 
    struct epsilon eps_c = cold_epsilon(env,K,w), eps_s = source_epsilon(env,K,w);
    double complex eps1 = 1. + eps_c.eps1 + eps_s.eps1, eps2 = eps_c.eps2 + eps_s.eps2;
    return ((eps1-n_pl*n_pl)/eps1 - n_pr*n_pr)*((eps1*eps1-eps2*eps2)/eps1 - n_pl*n_pl - n_pr*n_pr) - n_pl*n_pl*eps2/eps1*eps2/eps1;
    //return ((eps1*eps1-eps2*eps2)/eps1 - n_pr*n_pr);
}

static inline 
double N_dispersion_relation(double n, VectorH K, double const * w_vec) {
    double n_pl = K.pl/w_vec[0], n_pr = K.pr/w_vec[0];
    return n*n - n_pl*n_pl - n_pr*n_pr;
}


void cold_dispersion_relationH(void const *ptr, VectorSp R, VectorH K, double const * w_vec, double * res_vec) {
    epsilon_context_t const *cntx = ptr;
    PhysicalEnvironment phys_env = cntx->phys_env(cntx->phys_env_cntx,R);
    double complex w = w_vec[0] + I*w_vec[1];
    double complex res = cold_dispersion_relation(phys_env,K,w);
    res_vec[0] = creal(res); res_vec[1] = cimag(res);
}

void cold_dispersion_relationSp(void const *ptr, VectorSp R, VectorSp K, double const * w_vec, double * res_vec) {
    epsilon_context_t const *cntx = ptr;
    PhysicalEnvironment phys_env = cntx->phys_env(cntx->phys_env_cntx,R);
    double complex w = w_vec[0] + I*w_vec[1];
    double complex res = cold_dispersion_relation(phys_env,projection_of_on(K,phys_env.H),w);
    res_vec[0] = creal(res); res_vec[1] = cimag(res);
}


void source_dispersion_relationH(void const *ptr, VectorSp R, VectorH K, double const * w_vec, double * res_vec) {
    epsilon_context_t const *cntx = ptr;
    PhysicalEnvironment phys_env = cntx->phys_env(cntx->phys_env_cntx,R);
    double complex w = w_vec[0] + I*w_vec[1];
    double complex res = source_dispersion_relation(phys_env,K,w);
    res_vec[0] = creal(res); res_vec[1] = cimag(res);
}

void source_dispersion_relationSp(void const *ptr, VectorSp R, VectorSp K, double const * w_vec, double * res_vec) {
    epsilon_context_t const *cntx = ptr;
    PhysicalEnvironment phys_env = cntx->phys_env(cntx->phys_env_cntx,R);
    double complex w = w_vec[0] + I*w_vec[1];
    double complex res = source_dispersion_relation(phys_env,projection_of_on(K,phys_env.H),w);
    res_vec[0] = creal(res); res_vec[1] = cimag(res);
}


void warm_dispersion_relationH(void const *ptr, VectorSp R, VectorH K, double const * w_vec, double * res_vec) {
    epsilon_context_t const *cntx = ptr;
    PhysicalEnvironment phys_env = cntx->phys_env(cntx->phys_env_cntx,R);
    double complex w = w_vec[0] + I*w_vec[1];
    double complex res = warm_dispersion_relation(phys_env,K,w);
    res_vec[0] = creal(res); res_vec[1] = cimag(res);
}

void warm_dispersion_relationSp(void const *ptr, VectorSp R, VectorSp K, double const *w_vec, double *res_vec) {
    epsilon_context_t const *cntx = ptr;
    PhysicalEnvironment phys_env = cntx->phys_env(cntx->phys_env_cntx,R);
    double complex w = w_vec[0] + I*w_vec[1];
    double complex res = warm_dispersion_relation(phys_env,projection_of_on(K,phys_env.H),w);
    res_vec[0] = creal(res); res_vec[1] = cimag(res);
}


void N_dispersion_relationH(void const *ptr, VectorSp R, VectorH K, double const * w_vec, double * res_vec) {
    N_context_t const *cntx = ptr;
    res_vec[1] = 0.;
    res_vec[0] = N_dispersion_relation(cntx->n,K,w_vec);
}

void N_dispersion_relationSp(void const *ptr, VectorSp R, VectorSp K, double const * w_vec, double * res_vec) {
    N_context_t const *cntx = ptr;
    PhysicalEnvironment phys_env = cntx->phys_env(cntx->phys_env_cntx,R);

    res_vec[1] = 0.;
    res_vec[0] = N_dispersion_relation(cntx->n,projection_of_on(K,phys_env.H),w_vec);
}



