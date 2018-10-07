#include <broydn_solver.h>

#include <math.h>

/************
 * J[0]  J[1]
 * J[2]  J[3]
 ************/ 
void calculate_jacobian(func_2d F, void const * Fcntx, double const * x0, double dx, double dy, double * J) {
    double x_pd[2], x_md[2]; //x0 + dx -> pd  x0 - dx -> md
    double F_pd[2], F_md[2]; //F(x0+dx) and F(x0-dx)

    x_pd[0] = x0[0]+dx; x_pd[1] = x0[1];
    x_md[0] = x0[0]-dx; x_md[1] = x0[1];
    F(Fcntx,x_pd,F_pd); F(Fcntx,x_md,F_md);
    
    J[0] = 0.5*(F_pd[0]-F_md[0])/dx; J[2] = 0.5*(F_pd[1]-F_md[1])/dx;

    x_pd[0] = x0[0]; x_pd[1] = x0[1]+dy;
    x_md[0] = x0[0]; x_md[1] = x0[1]-dy;
    F(Fcntx,x_pd,F_pd); F(Fcntx,x_md,F_md);

    J[1] =  0.5*(F_pd[0]-F_md[0])/dy; J[3] = 0.5*(F_pd[1]-F_md[1])/dy; 
}

double det(double const * J) {
    return J[0]*J[3] - J[1]*J[2];
}

void invert_jacobian(double const * J, double detOver, double * inverted_J) {
    inverted_J[0] = J[3]*detOver; inverted_J[1] = - J[1]*detOver; inverted_J[2] = - J[2]*detOver; inverted_J[3] = J[0]*detOver;
}

void matrix_vector_product(double const * J, double const * vec, double * res_vec) {
    res_vec[0] = J[0]*vec[0] + J[1]*vec[1];
    res_vec[1] = J[2]*vec[0] + J[3]*vec[1];
}

void vector_tvector_product(double const *vec, double const *tvec, double * J) {
    J[0] = vec[0]*tvec[0]; J[1] = vec[0]*tvec[1];
    J[2] = vec[1]*tvec[0]; J[3] = vec[1]*tvec[1];
}

void matrix_solve(double const * J, double detOver, double const * x0_vec, double const * f0_vec, double * res_vec) {
    double inverted_J[4];
    invert_jacobian(J,detOver,inverted_J);

    double dx_vec[2];
    matrix_vector_product(inverted_J,f0_vec,dx_vec);

    res_vec[0] = x0_vec[0] - dx_vec[0];
    res_vec[1] = x0_vec[1] - dx_vec[1];
}

IterativeSolverStatus broydn_solver(
    func_2d F, void const * Fcntx, double const * initial_guess
    , double dx, double dy
    , double eps, unsigned int max_iteration
    , double * res_vec, unsigned int * res_iteration
) {
    IterativeSolverStatus res_status = SLV_OK;
    unsigned int iter_made = 0;
    
    double curr[2], Fcurr[2], next[2], Fnext[2];
    curr[0] = initial_guess[0]; curr[1] = initial_guess[1];
    next[0] = curr[0]; next[1] = curr[1];
    
    F(Fcntx,curr,Fcurr);
    Fnext[0] = Fcurr[0]; Fnext[1] = Fcurr[1];

    double J[4];
    calculate_jacobian(F,Fcntx,curr,dx,dy,J);

    while (hypot(Fnext[0],Fnext[1]) > eps) {
        
        if (iter_made > max_iteration) {
            res_status = SLV_MAX_ITER;
            break;
        }

        double detOver = 1.0/det(J);
        if( isnan(detOver) || isinf(detOver) ) {
            res_status = SLV_SINGULARITY;
            break;
        }
        
        matrix_solve(J,detOver,curr,Fcurr,next);
        F(Fcntx,next,Fnext);

        //recalculate Jacobian matrix
        //J(j+1) = J(j) + (dy - J*dx)dxT/(dx*dxT)
        //dx = x(j+1) - x(j); dy = F(j+1) - F(j)
        double dx_vec[2] = {next[0] - curr[0], next[1] - curr[1]};
        double dy_vec[2] = {Fnext[0] - Fcurr[0], Fnext[1] - Fcurr[1]};

        double tmp_vec[2];
        matrix_vector_product(J,dx_vec,tmp_vec);

        dy_vec[0] -= tmp_vec[0]; dy_vec[1] -= tmp_vec[1];

        double tmp_J[4];
        vector_tvector_product(dy_vec,dx_vec,tmp_J);

        double dx_norm = dx_vec[0]*dx_vec[0] + dx_vec[1]*dx_vec[1];
        tmp_J[0] /= dx_norm; tmp_J[1] /= dx_norm; tmp_J[2] /= dx_norm; tmp_J[3] /= dx_norm;

        J[0] += tmp_J[0]; J[1] += tmp_J[1]; J[2] += tmp_J[2]; J[3] += tmp_J[3];

        ++iter_made;
        curr[0] = next[0]; curr[1] = next[1];
        Fcurr[0] = Fnext[0]; Fcurr[1] = Fnext[1];
    }

    *res_iteration = iter_made;
    res_vec[0] = next[0];
    res_vec[1] = next[1];

    return res_status;
}


