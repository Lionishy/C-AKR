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
        , .cold_density = 0.999, .source_density = 0.001
    };

    FILE *map_file_up = fopen("./new_map_up1.txt","w");
    FILE *map_file_dn = fopen("./new_map_dn1.txt","w");
    double dw = 0.00001;

    for (double w = 2.0; w > 0.; w -= dw) {
        VectorH K = {0.,0.};
        double res = warm_left_dispersion_relation(
            homogeneous_physical_environment(&physical_environment_cntx,R0)
            , K
            , w
            );

        if (res > 0.)
            fprintf(map_file_up,"%f %f %f\n",w,sqrt(res),sqrt(res)*w);
        else
            printf("Negative %f %f\n",w,res);
    }

    /*for (double w = 1.+dw; w > 0.; w -= dw) {
        VectorH K = {0.,0.};
        double res = warm_left_dispersion_relation(
            homogeneous_physical_environment(&physical_environment_cntx,R0)
            , K
            , w
            );

        if (res > 0.)
            fprintf(map_file_dn,"%f %f %f\n",w,sqrt(res),sqrt(res)*w);
    }*/

    END: {
        if (NULL != map_file_dn) fclose(map_file_dn);
        if (NULL != map_file_up) fclose(map_file_up);
    }

}