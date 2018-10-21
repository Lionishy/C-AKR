#include <physical_environment.h>
#include <new_dispersion_relation.h>
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
    homogeneous_physical_environment_context_t physical_environment_cntx = {
          .R0 = R0, .V0 = V0
        , .omega_cc0 = 1., .omega_pc0 = 0.1
        , .cold_density = 0.01, .source_density = 0.99
    };

    FILE *map_file = fopen("./new_map.txt","w");

    for (double w = 2.0; w > -2.0; w -= 0.0005) {
        VectorH K = {0.,0.};
        double complex w_cplx = w;
        double complex res = left_dispersion_relation(
            homogeneous_physical_environment(&physical_environment_cntx,R0)
            , K
            , w_cplx
            );

        fprintf(map_file,"%f %f %f\n",w,creal(res),cimag(res));
    }

    END: {
        if (NULL != map_file) fclose(map_file);
    }

}