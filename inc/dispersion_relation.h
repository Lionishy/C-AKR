#ifndef IKI_DispersionRelation_H
#define IKI_DispersionRelation_H

#include <physical_environment.h>
#include <vector.h>

#include <complex.h>

typedef void (*relationSp_t) (void const *cntx, VectorSp R, VectorSp K, double const *w_vec, double *res_vec);
typedef void (*relationH_t) (void const *cntx, VectorSp R, VectorH K, double const *w_vec, double *res_vec);

typedef struct _epsilon_context {
    phys_env_t phys_env;
    void const *phys_env_cntx;
} epsilon_context_t;

void cold_dispersion_relationSp(void const *cntx, VectorSp R, VectorSp K, double const *w_vec, double *res_vec);
void cold_dispersion_relationH(void const *cntx, VectorSp R, VectorH K, double const *w_vec, double *res_vec);

void source_dispersion_relationSp(void const *cntx, VectorSp R, VectorSp K, double const *w_vec, double *res_vec);
void source_dispersion_relationH(void const *cntx, VectorSp R, VectorH K, double const *w_vec, double *res_vec);

void warm_dispersion_relationSp(void const *cntx, VectorSp R, VectorSp K, double const *w_vec, double *res_vec);
void warm_dispersion_relationH(void const *cntx, VectorSp R, VectorH K, double const *w_vec, double *res_vec);


typedef struct _N_context {
    phys_env_t phys_env;
    void const *phys_env_cntx;
    double n;
} N_context_t;

void N_dispersion_relationSp(void const *cntx, VectorSp R, VectorSp K, double const *w_vec, double *res_vec);
void N_dispersion_relationH(void const *cntx, VectorSp R, VectorH K, double const *w_vec, double *res_vec);

#endif //IKI_DispersionRelation_H