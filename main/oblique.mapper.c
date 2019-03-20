#include <physical_environment.h>
#include <dispersion_relation.h>
#include <secant_solver.h>
#include <vector.h>

#include <stdbool.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

typedef struct _oblique_real_K_context {
    double w, Kpl, Kpr;
    VectorSp R;
    void const *dispersion_cntx;
} oblique_real_K_cntx_t;

void oblique_real_K(void const *ptr, double K, double *out) {
    oblique_real_K_cntx_t const *cntx = ptr;
    double res[2]; double w_vec[2] = {cntx->w,0.};
    warm_dispersion_relation_minusH(cntx->dispersion_cntx,cntx->R,(VectorH){cntx->Kpl*K,cntx->Kpr*K},w_vec,res);
    *out = res[0];
}

enum WriteState {WAIT,WRITE};

int main() {
    //ФИЗИЧЕСКАЯ МОДЕЛЬ
    VectorSp R0 = {1.975,0.44928364,0.}; //точка нормировки
    VectorH V0  = {0.05,0.1};

    //параметры окружающей плазмы
    homogeneous_physical_environment_context_t physical_environment_cntx = {
          .R0 = R0, .V0 = V0
        , .omega_cc0 = 1., .omega_pc0 = 0.1
        , .cold_density = 0.064, .source_density = 1.0
    };

    //контекст для дисперсионного уравнения
    epsilon_context_t dispersion_relation_cntx = {
        homogeneous_physical_environment, &physical_environment_cntx
    };


    double angle = 3.14159265358979323846/3.;
    //контекст для действительной функции от модуля К
    oblique_real_K_cntx_t oblique_cntx = {
        2.0, -sin(angle), cos(angle)
        , R0
        , &dispersion_relation_cntx
    };    

    FILE * fd = fopen("./oblique.txt","w");
    if (NULL == fd) {
        printf("Can't open file!\n");
        return 0;
    }

    unsigned count = 0; double K1 = 0.001, K2 = 2.0, w1 = 0.001, w2 = 2.0;
    double K = K1;
    for (double w = w1; w < 2.0; w += 1.e-6) {
        double res = -1; unsigned iterations;
        
        oblique_cntx.w = w;
        IterativeSolverStatus status =  secant_solver(
            oblique_real_K, &oblique_cntx
            , K
            , 1.0e-7, 1.0e-9, 1000
            , &res, &iterations
        );

        if (SLV_OK == status) {
            K = res;
            if (0 == (count++%50) && K >= 0.) {
                fprintf(fd,"%f %f %f\n",w,K/w,K);     
                count = 1;
            }
            K *= K < 0 ? -1. : 1.;
        }
    } 

    fclose(fd);


    return 0;
}
