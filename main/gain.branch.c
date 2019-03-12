#include <physical_environment.h>
#include <dispersion_relation.h>
#include <correction.h>
#include <vector.h>

#include <stdbool.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

int main() {
    //ФИЗИЧЕСКАЯ МОДЕЛЬ
    VectorSp R0 = {1.975,0.44928364,0.}; //точка нормировки
    //VectorH V0  = {0.05,0.12247449};
    VectorH V0 = {0.05,0.1};

    //параметры окружающей плазмы
    /*dipole_physical_environment_context_t physical_environment_cntx = {
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

    //контекс для дисперсионного уравнения
    /*epsilon_context_t dispersion_relation_cntx = {
        dipole_physical_environment, &physical_environment_cntx
    };*/
    epsilon_context_t dispersion_relation_cntx = {
        homogeneous_physical_environment, &physical_environment_cntx
    };


    omega_correctorH_context_t corrector_cntx = {
          1.e-9, 1.e-9, 1.e-10, 1000u
        , warm_dispersion_relation_minusH, &dispersion_relation_cntx
    };

    //start point
    unsigned counter = 0;
    VectorSp R = R0;
    double angle = 3.14159265358979323846/6.;
    VectorH K = {sin(angle),cos(angle)}; double kmod = hypot(K.pl,K.pr); K.pl /= kmod; K.pr /= kmod; 
    double w_vec[2] = { 0.99876, 0.00498}; double n_start = 0; //0.99875607, 0.00497512

    FILE *fd = fopen("./gain.branch.txt","w");
    for (double n = n_start; n < 2.0; n += 1.e-6) {
        VectorH Kh = (VectorH){K.pl*n, K.pr*n};
        if (omega_correctorH(&corrector_cntx,&R0,&Kh,w_vec)) {
            w_vec[1] = w_vec[1] < 0. ? -w_vec[1] : w_vec[1];
            if (w_vec[0] < 0.) {
                printf("Negative omega\n");
                printf("%.8f %.8f %.8f\n",n,w_vec[0],w_vec[1]);
                return 0;
            }
            if (0 == (counter %= 100))fprintf(fd,"%.8f %.8f %.8f %.8f\n",w_vec[0],hypot(Kh.pl,Kh.pr)/w_vec[0],hypot(Kh.pl,Kh.pr),w_vec[1]); 
            ++counter;
        }
        else {
            printf("%f Error\n",n); fflush(stdout);
            return 0;
        }
    }

    fclose(fd);

    return 0;
}


