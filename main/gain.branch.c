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
    VectorH V0  = {0.05,0.12247449};

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
        , .cold_density = 0.01, .source_density = 0.99
    };

    //контекс для дисперсионного уравнения
    /*epsilon_context_t dispersion_relation_cntx = {
        dipole_physical_environment, &physical_environment_cntx
    };*/
    epsilon_context_t dispersion_relation_cntx = {
        homogeneous_physical_environment, &physical_environment_cntx
    };


    omega_correctorH_context_t corrector_cntx = {
          1.e-8, 1.e-8, 1.e-12, 1000u
        , warm_dispersion_relationH, &dispersion_relation_cntx
    };

    //start point
    unsigned counter = 0;
    VectorSp R = R0;
    VectorH K = {1.,1.}; double kmod = hypot(K.pl,K.pr); K.pl /= kmod; K.pr /= kmod; 
    double w_vec[2] = {1.03580441,0.}; double n_start = 1.2; // 1.03580441

    FILE *fd = fopen("./gain.branch.txt","w");
    for (double n = n_start; n > 0.1; n -= 1.e-7) {
        VectorH Kh = (VectorH){K.pl*n, K.pr*n};
        if (omega_correctorH(&corrector_cntx,&R0,&Kh,w_vec)) {
            w_vec[1] = w_vec[1] < 0. ? -w_vec[1] : w_vec[1];
            if (w_vec[0] < 0.) {
                printf("Negative omega\n");
                printf("%.8f %.8f %.8f\n",n,w_vec[0],w_vec[1]);
                return 0;
            }
            if (0 == (counter %= 50))fprintf(fd,"%.8f %.8f %.8f %.8f\n",n,w_vec[0],w_vec[1],hypot(Kh.pl,Kh.pr)/w_vec[0]); 
            ++counter;
        }
        else {
            printf("%f Error\n",n); fflush(stdout);
            return 0;
        }
    }

    /*for (double n = n_start; n < 1.0; n += 1.e-5) {
        VectorH Kh = (VectorH){K.pl*n, K.pr*n};
        if (omega_correctorH(&corrector_cntx,&R0,&Kh,w_vec)) {
            if (w_vec[0] < 0.) {
                printf("Negative omega\n");
                printf("%.8f %.8f %.8f\n",n,w_vec[0],w_vec[1]);
                return 0;
            }
            fprintf(fd,"%.8f %.8f %.8f %.8f\n",n,w_vec[0],w_vec[1],n/w_vec[0]); 
        }
        else {
            printf("Error\n");
            return 0;
        }
    }*/

    fclose(fd);

    return 0;
}


