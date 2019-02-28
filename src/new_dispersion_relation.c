#include <physical_environment.h>

#include <complex.h>
#include <math.h>

struct epsilon {
    double eps1, eps2;
};

static inline 
struct epsilon source_epsilon(PhysicalEnvironment env, VectorH K, double w) {
    double w_sqr = w*w;
    double omega_pss = env.omega_ps*env.omega_ps;

    double sf = 1. - 0.5*(env.V.pl*env.V.pl + env.V.pr*env.V.pr);
    double kv_pl = K.pl*env.V.pl;

    double
        f1 = w-env.omega_cc*sf-kv_pl,
        f2 = w+env.omega_cc*sf-kv_pl;

    double
        F1 = (w-kv_pl)/f1 + 0.5*env.V.pr*env.V.pr*(K.pl*K.pl-w*env.omega_cc)/(f1*f1),
        F2 = (w-kv_pl)/f2 + 0.5*env.V.pr*env.V.pr*(K.pl*K.pl-w*env.omega_cc)/(f2*f2);

    double
        eps1 = -0.5*omega_pss/w_sqr*(F1 + F2),
        eps2 =  0.5*omega_pss/w_sqr*(F1 - F2);

    return (struct epsilon){eps1,eps2};
}

static inline 
struct epsilon cold_epsilon(PhysicalEnvironment env, VectorH K, double w) {
    double w_sqr = w*w;
    double omega_pcs = env.omega_pc*env.omega_pc;
    double omega_ccs = env.omega_cc*env.omega_cc;

    double 
        eps1 = -(omega_pcs)/(w_sqr-omega_ccs),
        eps2 = env.omega_cc/w*omega_pcs/(w_sqr-omega_ccs);

   return (struct epsilon){eps1,eps2}; 
}

/**
 * Asking for a real frequency (omega)
 * And returning a squared refractive index (Npr*Npr)
*/
double warm_plus_Npr(PhysicalEnvironment env, VectorH K, double w) {
    double n_pl = K.pl/w, n_pr = K.pr/w; 
    struct epsilon eps_c = cold_epsilon(env,K,w), eps_s = source_epsilon(env,K,w);
    double eps1 = 1. + eps_c.eps1 + eps_s.eps1, eps2 = eps_c.eps2 + eps_s.eps2; 
    double alpha = (eps1 - n_pl*n_pl)/eps1;
    double betta = (eps1*eps1 - eps2*eps2)/eps1 - n_pl*n_pl;
    double gamma = alpha-betta;
    return (alpha+betta)*0.5 + sqrt(gamma*gamma + 4.0*(n_pl*eps2/eps1)*(n_pl*eps2/eps1))*0.5*(gamma > 0 ? 1. : -1.);
}

double warm_minus_Npr(PhysicalEnvironment env, VectorH K, double w) {
    double n_pl = K.pl/w, n_pr = K.pr/w; 
    struct epsilon eps_c = cold_epsilon(env,K,w), eps_s = source_epsilon(env,K,w);
    double eps1 = 1. + eps_c.eps1 + eps_s.eps1, eps2 = eps_c.eps2 + eps_s.eps2; 
    double alpha = (eps1 - n_pl*n_pl)/eps1;
    double betta = (eps1*eps1 - eps2*eps2)/eps1 - n_pl*n_pl;
    double gamma = alpha-betta;
    return (alpha+betta)*0.5 + sqrt((alpha-betta)*(alpha-betta) + 4.0*(n_pl*eps2/eps1)*(n_pl*eps2/eps1))*0.5*(gamma > 0 ? -1. : 1.);
}

