#include <physical_environment.h>

#include <complex.h>
#include <math.h>

struct epsilon {
    double complex eps1, eps2, f0;
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
        F1 = (w-kv_pl)/f1 + 0.5*env.V.pr*env.V.pr*(K.pl*K.pl-w*env.omega_cc)/(f1*f1),
        F2 = (w-kv_pl)/f2 + 0.5*env.V.pr*env.V.pr*(K.pl*K.pl-w*env.omega_cc)/(f2*f2);

    double complex
        eps1 = -0.5*omega_pss/w_sqr*(F1 + F2),
        eps2 =  0.5*omega_pss/w_sqr*(F1 - F2);

    return (struct epsilon){eps1,eps2,f2*f2*f1*f1};
}

static inline 
struct epsilon cold_epsilon(PhysicalEnvironment env, VectorH K, double complex w) {
    double complex w_sqr = w*w;
    double omega_pcs = env.omega_pc*env.omega_pc;
    double omega_ccs = env.omega_cc*env.omega_cc;

    double f = (w_sqr-omega_ccs);

    double complex
        eps1 = -(omega_pcs)/f,
        eps2 = env.omega_cc/w*omega_pcs/f;

   return (struct epsilon){eps1,eps2,f}; 
}

double complex left_dispersion_relation(PhysicalEnvironment env, VectorH K, double complex w) {
    struct epsilon eps_c = cold_epsilon(env,K,w), eps_s = source_epsilon(env,K,w);
    double complex 
        eps1 = 1. + eps_c.eps1 + eps_s.eps1,
        eps2 = eps_c.eps2 + eps_s.eps2;

    return sqrt((eps1*eps1-eps2*eps2)/eps1);
}
