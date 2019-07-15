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

typedef struct _real_Kp_context {
    double w, Kpr;
    VectorSp R;
    void const *dispersion_cntx;
} real_Kpl_cntx_t;

void real_Kpl(void const *ptr, double Kpl, double *out) {
    real_Kpl_cntx_t const *cntx = ptr;
    double res[2]; double w_vec[2] = {cntx->w,0.};
    warm_dispersion_relation_minusH(cntx->dispersion_cntx,cntx->R,(VectorH){Kpl,cntx->Kpr},w_vec,res);
    *out = res[0];
}

enum WriteState {WAIT,WRITE};

int main() {
    //ФИЗИЧЕСКАЯ МОДЕЛЬ
    VectorSp R0 = {1.975,0.44928364,0.}; //точка нормировки
    VectorH V0  = {0.05,0.15};

    //параметры окружающей плазмы
    homogeneous_physical_environment_context_t physical_environment_cntx = {
          .R0 = R0, .V0 = V0
        , .omega_cc0 = 1., .omega_pc0 = 0.1
        , .cold_density = 0.001, .source_density = 0.2
    };

    //контекст для дисперсионного уравнения
    epsilon_context_t dispersion_relation_cntx = {
        homogeneous_physical_environment, &physical_environment_cntx
    };

    //контекст для действительной функции от модуля К
    real_Kpl_cntx_t Kpl_cntx = {
        1.0, 0.8
        , R0
        , &dispersion_relation_cntx
    };    

    FILE * fd = fopen("./kpl.map.txt","w");
    if (NULL == fd) {
        printf("Can't open file!\n");
        return 0;
    }

    unsigned count = 0; double Kpl1 = 0.001, Kpl2 = 2.0, w1 = Kpl_cntx.Kpr, w2 = 2.0;
    double Kpl = Kpl1;
    for (double w = w1; w < 2.0; w += 1.e-6) {
        double res = -1; unsigned iterations;
        
        Kpl_cntx.w = w;
        IterativeSolverStatus status =  secant_solver(
            real_Kpl, &Kpl_cntx
            , Kpl
            , 1.0e-7, 1.0e-9, 1000
            , &res, &iterations
        );

        if (SLV_OK == status) {
            Kpl = res;
            if (0 == (count++%50) && Kpl >= 0.) {
                fprintf(fd,"%f %f %f\n",w,hypot(Kpl,Kpl_cntx.Kpr)/w,Kpl);     
                count = 1;
            }
            Kpl *= Kpl < 0 ? -1. : 1.;
        }
    } 

    fclose(fd);


    return 0;
}
