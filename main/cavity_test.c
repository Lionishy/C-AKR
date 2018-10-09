#include <vector.h>

#include <stdio.h>
#include <math.h>


typedef struct _tanhcavity_cntx {
    double L0, widthL, dL;
} tanhcavity_cntx_t;

double tanhcavity(void *ptr, double L) {
    tanhcavity_cntx_t *cntx = ptr;
    double rightL = cntx->L0+cntx->widthL, leftL = cntx->L0-cntx->widthL;

    return 1. + 0.5*(tanh((L-rightL)/cntx->dL) + tanh((leftL-L)/cntx->dL));
}

typedef struct _expcavity_cntx {
    double L0, widthL;
} expcavity_cntx_t;

double expcavity(void *ptr, double L) {
    expcavity_cntx_t *cntx = ptr;
    return exp( -(L-cntx->L0)*(L-cntx->L0)/(cntx->widthL*cntx->widthL) );
}


int main() {
    tanhcavity_cntx_t tanhcavity_cntx = {
        10.45, 0.01, 0.001
    };

    expcavity_cntx_t expcavity_cntx = {
        10.45, 0.005
    };
    
    FILE *cav = fopen("./cavity-test.txt","w");
    if (NULL == cav) {
        printf("Can't open file ../cavity-test.txt to write result\n");
        goto END;
    }

    for (double start = tanhcavity_cntx.L0-2*tanhcavity_cntx.widthL, end = tanhcavity_cntx.L0+2*tanhcavity_cntx.widthL; start < end; start += tanhcavity_cntx.widthL/100.) {
        double thcavity_res = tanhcavity(&tanhcavity_cntx,start);
        double excavity_res = 1. - expcavity(&expcavity_cntx,start);
        fprintf(cav,"%f %.9f %.9f\n",start,thcavity_res,excavity_res);
    } 

    END:
    if (NULL != cav) fclose(cav);

    return 0;
}