#include <physical_environment.h>
#include <dispersion_relation.h>
#include <secant_solver.h>
#include <vector.h>

#include <stdbool.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

typedef struct _real_Kr_context {
    double w, Kpl;
    VectorSp R;
    void const *dispersion_cntx;
} real_Kpr_cntx_t;

void real_Kpr(void const *ptr, double Kpr, double *out) {
    real_Kpr_cntx_t const *cntx = ptr;
    double res[2]; double w_vec[2] = {cntx->w,0.};
    warm_dispersion_relation_minusH(cntx->dispersion_cntx,cntx->R,(VectorH){cntx->Kpl,Kpr},w_vec,res);
    *out = res[0];
}

enum WriteState {WAIT,WRITE};

int main() {
    //ФИЗИЧЕСКАЯ МОДЕЛЬ
    VectorSp R0 = {1.975,0.44928364,0.}; //точка нормировки
    VectorH V0  = {0.0,0.15};

    //параметры окружающей плазмы
    homogeneous_physical_environment_context_t physical_environment_cntx = {
          .R0 = R0, .V0 = V0
        , .omega_cc0 = 1., .omega_pc0 = 0.1
        , .cold_density = 0.1, .source_density = 1.0
    };

    //контекст для дисперсионного уравнения
    epsilon_context_t dispersion_relation_cntx = {
        homogeneous_physical_environment, &physical_environment_cntx
    };

    //контекст для действительной функции от модуля К
    real_Kpr_cntx_t Kpr_cntx = {
        1.0, 0.
        , R0
        , &dispersion_relation_cntx
    };    

    FILE * fd = fopen("./kpr.map.txt","w");
    if (NULL == fd) {
        printf("Can't open file!\n");
        return 0;
    }

    unsigned count = 0; double Kpr1 = 0.001, Kpr2 = 2.0, w1 = Kpr_cntx.Kpl, w2 = 2.0;
    double Kpr = Kpr2;
    for (double w = w2; w > 0.0; w -= 1.e-6) {
        double res = -1; unsigned iterations;
        
        Kpr_cntx.w = w;
        IterativeSolverStatus status =  secant_solver(
            real_Kpr, &Kpr_cntx
            , Kpr
            , 1.0e-7, 1.0e-9, 1000
            , &res, &iterations
        );

        if (SLV_OK == status) {
            Kpr = res;
            if (0 == (count++%50) && Kpr >= 0.) {
                fprintf(fd,"%f %f %f\n",w,hypot(Kpr,Kpr_cntx.Kpl)/w,Kpr);     
                count = 1;
            }
            Kpr *= Kpr < 0 ? -1. : 1.;
        }
    } 

    fclose(fd);


    return 0;
}