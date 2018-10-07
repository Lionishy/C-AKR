#ifndef IKI_ray_trace_H
#define IKI_ray_trace_H

#include <dispersion_relation.h>
#include <correction.h>
#include <vector.h>

#include <stdio.h>

void dD_dK(relationSp_t relation, void const * relation_cntx, VectorSp R, VectorSp K, double const * w_vec, VectorSp dK, VectorSp * res);
void dK_dR(relationSp_t relation, void const * relation_cntx, VectorSp R, VectorSp K, double const * w_vec, VectorSp dR, VectorSp * res);
void velocityK(VectorSp R, VectorSp dD_dR, double dD_dw, VectorSp *VK); 
void velocityR(VectorSp R, VectorSp dD_dK, double dD_dw, VectorSp *VR);


typedef enum _step_status {STP_OK, STP_ERROR_CORRECTION, STP_ERROR_ANY} step_status_t;
typedef step_status_t (*step_t)(void const *cntx, VectorSp *R, VectorSp *K, double *w_vec, VectorSp *VR, VectorSp *VK, double dt);


typedef struct _predictor_step_context {
    VectorSp dR, dK;
    double dw;
    relationSp_t relation;
    void const *relation_cntx;
    correctorSp_t corrector;
    void const *corrector_cntx;
} predictor_step_context_t;

step_status_t predictor_step(void const *env_ptr, VectorSp *R, VectorSp *K, double *w_vec, VectorSp *VR, VectorSp *VK, double dt);


typedef struct _N1_step_env {
    VectorSp dR, dK;
    double dw;
    relationSp_t relation;
    void const *relation_cntx;
} N1_step_env_t;

step_status_t N1_step(void const *env_ptr, VectorSp * R, VectorSp * K, double * w_vec, VectorSp *VR, VectorSp *VK, double dt);


#endif