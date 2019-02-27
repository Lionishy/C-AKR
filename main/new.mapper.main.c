#include <physical_environment.h>
#include <new_dispersion_relation.h>
#include <vector.h>

#include <stdbool.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

enum WriteState {WAIT,WRITE};

int main() {
    //ФИЗИЧЕСКАЯ МОДЕЛЬ
    VectorSp R0 = {1.975,0.44928364,0.}; //точка нормировки
    VectorH V0  = {0.0,0.05};

    //параметры окружающей плазмы
    homogeneous_physical_environment_context_t physical_environment_cntx = {
          .R0 = R0, .V0 = V0
        , .omega_cc0 = 1., .omega_pc0 = 0.1
        , .cold_density = 0.001, .source_density = 1.0
    };


    double prevSolve[2] = {0.,0.}; 
    char name[] = "./new_map_branch_";
    char ext[] = ".txt";
    int branch_counter = 0;
    enum WriteState st = WAIT; 

    FILE *map_file = NULL;
    double dw = 0.00001;

    for (double w = 2.0; w > 0.; w -= dw) {
        VectorH K = {0.,0.};
        double res = warm_minus_Npr(
            homogeneous_physical_environment(&physical_environment_cntx,R0)
            , K
            , w
            );

        if (res > 0.) {
            if (WAIT == st) {
                char count_name[3];
                sprintf(count_name,"%d",++branch_counter);
                char* full_name = malloc(strlen(name)+strlen(count_name)+strlen(ext)+1);
                strcpy(full_name,name);
                strcat(full_name,count_name);
                strcat(full_name,ext);
                printf("new branch: %s\n",full_name);

                map_file = fopen(full_name,"w");
                free(full_name);
                if (NULL == map_file) {
                    printf("Can't open file...");
                    return 0;
                }

                if (res < 1. && prevSolve[1] < 0.) {
                    double wMid = prevSolve[0] - prevSolve[1]*(w-prevSolve[0])/(res-prevSolve[1]);
                    fprintf(map_file,"%f %f %f\n",wMid,0,0);
                }

                st = WRITE;
            }
            fprintf(map_file,"%f %f %f\n",w,sqrt(res),sqrt(res)*w);
        } else {
            if (WRITE == st) {
                if (prevSolve[1] < 1.) {
                    double wMid = prevSolve[0] - prevSolve[1]*(w-prevSolve[0])/(res-prevSolve[1]);
                    fprintf(map_file,"%f %f %f\n",wMid,0,0);
                }                

                fclose(map_file);
                st = WAIT;
            }
        }
        prevSolve[0] = w; prevSolve[1] = res;   
    }

    END: {
        if (NULL != map_file && WRITE == st) fclose(map_file);
    }

}