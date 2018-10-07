#ifndef IKI_Vector_H
#define IKI_Vector_H

typedef struct _VectorSp {
    double r,th,phi;
} VectorSp;

typedef struct _VectorH {
    double pl,pr;
} VectorH;

double vecsp_norm(VectorSp vec);
double vecsp_scalar_prod(VectorSp vec1, VectorSp vec2);
VectorH projection_of_on(VectorSp of, VectorSp on);
void projection_of_on_insp(VectorSp of, VectorSp on, VectorSp *pl, VectorSp *pr);
VectorSp rotate(VectorSp to_rotate, VectorSp R, VectorSp dR);

#endif //IKI_Vector_H