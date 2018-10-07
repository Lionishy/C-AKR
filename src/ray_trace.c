#include <ray_trace.h>
#include <correction.h>
#include <vector.h>

#include <stdbool.h>
#include <math.h>

    #include <stdio.h>

void dD_dw(relationSp_t relation, void const *relation_cntx, VectorSp R0, VectorSp K0, double const *w_vec, double dw, double *res) {
    double pw[2], mw[2]; //w+dw -> pw   w-dw -> mw
    double pRel[2], mRel[2]; //relation(w+dw) -> pRel   relation(w-dw) -> mRel

    pw[0] = w_vec[0]+dw; pw[1] = w_vec[1];
    mw[0] = w_vec[0]-dw; mw[1] = w_vec[1];

    relation(relation_cntx,R0,K0,pw,pRel);
    relation(relation_cntx,R0,K0,mw,mRel);

    *res = 0.5*(pRel[0]-mRel[0])/dw;
}

void dD_dK(relationSp_t relation, void const *relation_cntx, VectorSp R0, VectorSp K0, double const *w_vec, VectorSp dK, VectorSp *res) {
    VectorSp pdK, mdK; //K+dK -> pdK   K-dK -> mdK
    double pRel[2], mRel[2]; //relation(K+dK) -> pRel   relation(K-dK) -> mRel
    
    //r
    pdK.r = K0.r+dK.r; pdK.th = K0.th; pdK.phi = K0.phi;
    mdK.r = K0.r-dK.r; mdK.th = K0.th; mdK.phi = K0.phi;
    relation(relation_cntx, R0, pdK, w_vec, pRel);
    relation(relation_cntx, R0, mdK, w_vec, mRel);
    res->r = 0.5*(pRel[0] - mRel[0])/dK.r;

    //th
    pdK.r = K0.r; pdK.th = K0.th+dK.th; pdK.phi = K0.phi;
    mdK.r = K0.r; mdK.th = K0.th-dK.th; mdK.phi = K0.phi;
    relation(relation_cntx, R0, pdK, w_vec, pRel);
    relation(relation_cntx, R0, mdK, w_vec, mRel);
    res->th = 0.5*(pRel[0] - mRel[0])/dK.th;

    //phi
    pdK.r = K0.r; pdK.th = K0.th; pdK.phi = K0.phi+dK.phi;
    mdK.r = K0.r; mdK.th = K0.th; mdK.phi = K0.phi-dK.phi;
    relation(relation_cntx, R0, pdK, w_vec, pRel);
    relation(relation_cntx, R0, mdK, w_vec, mRel);
    res->phi = 0.5*(pRel[0] - mRel[0])/dK.phi;
}

void dD_dR(relationSp_t relation, void const *relation_cntx, VectorSp R0, VectorSp K0, double const *w_vec, VectorSp dR, VectorSp *res) {
    VectorSp pdR, mdR; //R+dR -> pdR   R-dR -> mdR
    double pRel[2], mRel[2]; //relation(R+dR) -> pRel   relation(R-dR) -> mRel
    
    //r
    pdR.r = R0.r+dR.r; pdR.th = R0.th; pdR.phi = R0.phi;
    mdR.r = R0.r-dR.r; mdR.th = R0.th; mdR.phi = R0.phi;
    relation(relation_cntx, pdR, K0, w_vec, pRel);
    relation(relation_cntx, mdR, K0, w_vec, mRel);
    res->r = 0.5*(pRel[0] - mRel[0])/dR.r;

    //th
    pdR.r = R0.r; pdR.th = R0.th+dR.th; pdR.phi = R0.phi;
    mdR.r = R0.r; mdR.th = R0.th-dR.th; mdR.phi = R0.phi;
    relation(relation_cntx, pdR, K0, w_vec, pRel);
    relation(relation_cntx, mdR, K0, w_vec, mRel);
    res->th = 0.5*(pRel[0] - mRel[0])/dR.th;

    //phi
    pdR.r = R0.r; pdR.th = R0.th; pdR.phi = R0.phi+dR.phi;
    mdR.r = R0.r; mdR.th = R0.th; mdR.phi = R0.phi-dR.phi;
    relation(relation_cntx, pdR, K0, w_vec, pRel);
    relation(relation_cntx, mdR, K0, w_vec, mRel);
    res->phi = 0.5*(pRel[0] - mRel[0])/dR.phi;
}

double dD_dw_correction(double dD_dw_val, double cutoff) {
    if (dD_dw_val > 0 && dD_dw_val < cutoff) dD_dw_val = cutoff;
    if (dD_dw_val < 0 && dD_dw_val > -cutoff) dD_dw_val = -cutoff;
    return dD_dw_val;
}

void velocityK(VectorSp R, VectorSp dD_dR, double dD_dw, VectorSp *VK) {
    VK->r = dD_dR.r/dD_dw;
    VK->th = dD_dR.th/dD_dw/R.r;
    VK->phi = dD_dR.phi/dD_dw/R.r/sin(R.th);
}

void velocityR(VectorSp R, VectorSp dD_dK, double dD_dw, VectorSp *VR) {
    VR->r = -dD_dK.r/dD_dw;
    VR->th = -dD_dK.th/dD_dw/R.r;
    VR->phi = -dD_dK.phi/dD_dw/R.r/sin(R.th);
}

void VR_back_correction(VectorSp *V) {
    static bool init = true;
    static VectorSp V_prev = {0.,0.,0.}; static VectorSp coeff = {0.06,0.06,0.06};

    if (!init) {
       if (V_prev.r*V->r < 0.  && fabs(V_prev.r-V->r) > coeff.r) {
           V->r *= -1.;
       }

       if (V_prev.th*V->th < 0. && fabs(V_prev.th-V->th) > coeff.th) {
           V->th *= -1.;
       }

       if (V_prev.phi*V->phi < 0. && fabs(V_prev.phi-V->phi) > coeff.phi) {
           V->phi *= -1.;
       }

       V_prev = *V; 
       return;
    }

    V_prev = *V;
    init = false;
}

void VK_back_correction(VectorSp *V) {
    static bool init = true;
    static VectorSp V_prev = {0.,0.,0.}; static VectorSp coeff = {0.06,0.06,0.06};

    if (!init) {
       if (V_prev.r*V->r < 0.  && fabs(V_prev.r-V->r) > coeff.r) {
           V->r *= -1.;
       }

       if (V_prev.th*V->th < 0. && fabs(V_prev.th-V->th) > coeff.th) {
           V->th *= -1.;
       }

       if (V_prev.phi*V->phi < 0. && fabs(V_prev.phi-V->phi) > coeff.phi) {
           V->phi *= -1.;
       }

       V_prev = *V; 
       return;
    }

    V_prev = *V;
    init = false;
}

step_status_t predictor_step(void const *ptr, VectorSp *R, VectorSp *K, double *w_vec, VectorSp *VR, VectorSp *VK, double dt) {
    predictor_step_context_t const *cntx = ptr;

    /***********
     * calculate new K
     * *********/
    {
        VectorSp K_predicted, VK_predicted;
        double dD_dw_res;
        VectorSp dD_dR_res;
        
        //prediction
        dD_dR(cntx->relation,cntx->relation_cntx,*R,*K,w_vec,cntx->dR,&dD_dR_res);
        dD_dw(cntx->relation,cntx->relation_cntx,*R,*K,w_vec,cntx->dw,&dD_dw_res);

        velocityK(*R,dD_dR_res,dD_dw_res,VK);
        K_predicted.r = K->r + VK->r*dt;
        K_predicted.th = K->th + VK->th*dt;
        K_predicted.phi = K->phi + VK->phi*dt;

        //new step
        dD_dR(cntx->relation,cntx->relation_cntx,*R,K_predicted,w_vec,cntx->dR,&dD_dR_res);
        dD_dw(cntx->relation,cntx->relation_cntx,*R,K_predicted,w_vec,cntx->dw,&dD_dw_res);

        velocityK(*R,dD_dR_res,dD_dw_res,&VK_predicted);
        VK->r = 0.5*(VK->r + VK_predicted.r); VK->th = 0.5*(VK->th + VK_predicted.th); VK->phi = 0.5*(VK->phi + VK_predicted.phi);
        VK_back_correction(VK);

        VectorSp K_new = {K->r + VK->r*dt, K->th + VK->th*dt, K->phi + VK->phi*dt};

        if ( !(cntx->corrector(cntx->corrector_cntx,R,&K_new,w_vec)) ) {
            //printf("Kstep Error: %f %f %f %f %f %f %f\n",K_new.r,K_new.th,K_new.phi,VK->r,VK->th,VK->phi,dD_dw_res);
            return STP_ERROR_CORRECTION;
        }
        *K = K_new;
    }

    /*********
     * calculate new R
     * *******/
    {
        VectorSp R_predicted, R_new, K_new, VR_predicted;
        double dD_dw_res;
        VectorSp dD_dK_res;
        double w_vec_new[2];

        //prediction
        dD_dK(cntx->relation,cntx->relation_cntx,*R,*K,w_vec,cntx->dK,&dD_dK_res);
        dD_dw(cntx->relation,cntx->relation_cntx,*R,*K,w_vec,cntx->dw,&dD_dw_res);

        velocityR(*R,dD_dK_res,dD_dw_res,VR);
        R_predicted.r = R->r + VR->r*dt;
        R_predicted.th = R->th + VR->th*dt;
        R_predicted.phi = R->phi + VR->phi*dt;

        //rotate and correct
        K_new = rotate(*K,*R,(VectorSp){VR->r*dt,VR->th*dt,VR->phi*dt});
        w_vec_new[0] = w_vec[0]; w_vec_new[1] = w_vec[1]; 
        if ( !(cntx->corrector(cntx->corrector_cntx,&R_predicted,&K_new,w_vec_new)) ) {
            //printf("Rpredictor Error: %f %f %f %f %f %f %f\n",R_predicted.r, R_predicted.th, R_predicted.phi ,VR->r,VR->th,VR->phi,dD_dw_res);
            return STP_ERROR_CORRECTION;
        }

        //new step
        dD_dK(cntx->relation,cntx->relation_cntx,R_predicted,K_new,w_vec_new,cntx->dK,&dD_dK_res);
        dD_dw(cntx->relation,cntx->relation_cntx,R_predicted,K_new,w_vec_new,cntx->dw,&dD_dw_res);

        velocityR(R_predicted,dD_dK_res,dD_dw_res,&VR_predicted);
        VR->r = 0.5*(VR->r + VR_predicted.r); VR->th = 0.5*(VR->th + VR_predicted.th); VR->phi = 0.5*(VR->phi + VR_predicted.phi);
        VR_back_correction(VR);

        //rotation
        K_new = rotate(*K,*R,(VectorSp){VR->r*dt, VR->th*dt, VR->phi*dt});
        R->r += VR->r*dt; R->th += VR->th*dt; R->phi += VR->phi*dt;

        //correction
        if ( !(cntx->corrector(cntx->corrector_cntx,R,&K_new,w_vec)) ) {
            //printf("R result Error: %f %f %f %f %f %f %f\n",R->r, R->th, R->phi ,VR->r,VR->th,VR->phi,dD_dw_res);
            return STP_ERROR_CORRECTION;
        }
        *K = K_new;
    }

    return STP_OK;
}