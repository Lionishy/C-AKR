#ifndef IKI_correction_H
#define IKI_correction_H

#include <physical_environment.h>
#include <dispersion_relation.h>
#include <vector.h>

#include <stdbool.h>

typedef bool (*correctorSp_t)(void const *cntx, VectorSp *R, VectorSp *K, double *w_vec);
typedef bool (*correctorH_t)(void const *cntx, VectorSp *R, VectorH *K, double *w_vec);

bool trivial_corrector(void const *cntx, VectorSp *R, VectorSp *K, double *w_vec);


typedef struct _omega_correctorH_context {
    double dx, dy, eps;
    unsigned max_iter;
    relationH_t relation;
    void const *relation_cntx;
} omega_correctorH_context_t;

bool omega_correctorH(void const *cntx, VectorSp *R, VectorH *K, double *w_vec);

typedef struct _omega_correctorSp_context {
    double dx, dy, eps;
    unsigned max_iter;
    relationSp_t relation;
    void const *relation_cntx;
} omega_correctorSp_context_t;

bool omega_correctorSp(void const *cntx, VectorSp *R, VectorSp *K, double *w_vec);


typedef struct _gamma_kpr_correctorH_context {
    double dx, dy, eps;
    unsigned max_iter;
    relationH_t relation;
    void const *relation_cntx;
} gamma_kpr_correctorH_context_t;

bool gamma_kpr_correctorH(void const *cntx, VectorSp *R, VectorH *K, double *w_vec);

typedef struct _gamma_kpr_correctorSp_context {
    phys_env_t phys_env;
    void const *phys_env_cntx;
    double dx, dy, eps;
    unsigned max_iter;
    relationH_t relation;
    void const *relation_cntx;
} gamma_kpr_correctorSp_context_t;

bool gamma_kpr_correctorSp(void const *cntx, VectorSp *R, VectorSp *K, double *w_vec);


typedef struct _gamma_kpl_correctorH_context {
    double dx, dy, eps;
    unsigned max_iter;
    relationH_t relation;
    void const *relation_cntx;
} gamma_kpl_correctorH_context_t;

bool gamma_kpl_correctorH(void const *cntx, VectorSp *R, VectorH *K, double *w_vec);

typedef struct _gamma_kpl_correctorSp_context {
    phys_env_t phys_env;
    void const *phys_env_cntx;
    double dx, dy, eps;
    unsigned max_iter;
    relationH_t relation;
    void const *relation_cntx;
} gamma_kpl_correctorSp_context_t;

bool gamma_kpl_correctorSp(void const *cntx, VectorSp *R, VectorSp *K, double *w_vec);


typedef struct _kpl_correctorH_context {
    double dx, eps;
    unsigned max_iter;
    relationH_t relation;
    void const *relation_cntx;
} kpl_correctorH_context_t;

bool kpl_correctorH(void const *cntx, VectorSp *R, VectorH *K, double *w_vec);

typedef struct _kpl_correctorSp_context {
    phys_env_t phys_env;
    void const *phys_env_cntx;
    double dx, eps;
    unsigned max_iter;
    relationH_t relation;
    void const *relation_cntx;
} kpl_correctorSp_context_t;

bool kpl_correctorSp(void const *cntx, VectorSp *R, VectorSp *K, double *w_vec);

#endif //IKI_K_correction_H