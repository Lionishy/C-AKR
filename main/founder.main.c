#include <physical_environment.h>
#include <dispersion_relation.h>
#include <correction.h>

#include <stdbool.h>
#include <stdio.h>
#include <math.h>

struct w_gamma_domain {
    double w_start, w_stop, gamma_start, gamma_stop;
    double w_step, gamma_step;
};

int main() {
    /* Initial conditions */
    VectorSp R0 = {1.975,0.44928364,0.};
    //VectorH V0  = {0.05,0.12247449};
    VectorH V0 = {0.0,0.1};

    /*dipole_physical_environment_context_t physical_env_cntx = {
          .R0 = R0, .V0 = V0
        , .omega_cc0 = 1.0, .omega_pc0 = 0.1
        , .phi_width = 0.23, .L_width = 0.005
        , .source_density_coeff = 0.2 
    };*/
    
    homogeneous_physical_environment_context_t physical_environment_cntx = {
          .R0 = R0, .V0 = V0
        , .omega_cc0 = 1., .omega_pc0 = 0.1
        , .cold_density = 0., .source_density = 1.
    };

    /*epsilon_context_t dispersion_relation_cntx = {
        dipole_physical_environment, &physical_env_cntx
    };*/
    epsilon_context_t dispersion_relation_cntx = {
        homogeneous_physical_environment, &physical_environment_cntx
    };
    
    //w_start <--> w_stop 
    /*****************  //gamma_stop
     *               *
     *               *
     *               *
     *               *  //gamma_start
     * ***************/
    struct w_gamma_domain domain = {0.99,1.01,0.,0.006,1.e-4,1.e-4};
    omega_correctorH_context_t omega_correctorH_cntx = {
          1.e-9,1.e-9,1.e-12,1000u
        , warm_dispersion_relation_minusH, &dispersion_relation_cntx
    };

    VectorSp R = R0;
    double angle = 3.14159265358979323846/12.;
    VectorH K = {sin(angle),cos(angle)}; double n = 1.8;
    double mod = hypot(K.pl,K.pr); K.pl /= mod; K.pr /= mod;
    K.pr *= n; K.pl *= n;

    bool found = false;
    for (double w=domain.w_start; w < domain.w_stop && !found; w += domain.w_step) {
        for (double g=domain.gamma_start; g < domain.gamma_stop && !found; g += domain.gamma_step) {
            double w_vec[2] = {w,g};
            if( omega_correctorH(&omega_correctorH_cntx,&R,&K,w_vec) && fabs(w_vec[1]) > 1.e-6) {
                printf("%.8f %.8f %.8f, %.8f\n",w,g,w_vec[0],w_vec[1]); fflush(stdout);
                found = true;
            }    
        }
    }

    if (!found) {
        printf("Can't find root\n");
    }

    return 0;
}