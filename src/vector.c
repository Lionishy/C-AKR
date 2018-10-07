#include <vector.h>

#include <math.h>

double vecsp_norm(VectorSp vec) {
    return vec.r*vec.r + vec.th*vec.th + vec.phi*vec.phi;
}

double vecsp_scalar_prod(VectorSp vec1, VectorSp vec2) {
    return vec1.r*vec2.r + vec1.th*vec2.th + vec1.phi*vec2.phi;
}

VectorH projection_of_on(VectorSp of, VectorSp on) {
    double of_norm = vecsp_norm(of), on_norm = vecsp_norm(on);
    
    VectorH res;
    res.pl = vecsp_scalar_prod(of,on)/sqrt(on_norm);
    res.pr = sqrt(of_norm - res.pl*res.pl);
    return res;
}

void projection_of_on_insp(VectorSp of, VectorSp on, VectorSp *pl, VectorSp *pr) {
    double pl_coeff = vecsp_scalar_prod(of,on)/vecsp_norm(on);
    *pl = (VectorSp){pl_coeff*on.r,pl_coeff*on.th,pl_coeff*on.phi};
    *pr = (VectorSp){of.r - pl->r, of.th - pl->th, of.phi - pl->phi};
}

VectorSp rotate(VectorSp t, VectorSp R, VectorSp dR) {
    VectorSp tmp;
    tmp.r   = t.r*cos(dR.th)   +  t.th*sin(dR.th) + t.phi*sin(dR.phi)*sin(R.th);
    tmp.th  = t.th*cos(dR.th)  -  t.r*sin(dR.th) + t.phi*sin(dR.phi)*cos(R.th);
    tmp.phi = t.phi*cos(dR.phi) -  sin(dR.phi)*(t.r*sin(R.th) + t.th*cos(R.th));
    return tmp;
}

