#include <physical_environment.h>
#include <dispersion_relation.h>
#include <correction.h>
#include <vector.h>

#include <stdbool.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

struct map_domain {
    double kpl_start, kpl_stop, w_start, w_stop, kpl_step, w_step;
};

static inline
bool linear_search(double x, double y, double * points, unsigned count_points, double tolerance) {
    for (unsigned count = 0; count != count_points; ++count)
        if ( fabs( (x-points[2*count])/points[2*count] ) < tolerance && fabs( (y-points[2*count+1])/points[2*count+1] ) < tolerance  )
            return true;
    return false;
}

int main() {
    //ФИЗИЧЕСКАЯ МОДЕЛЬ
    VectorSp R0 = {1.975,0.44928364,0.}; //точка нормировки
    VectorH V0  = {0.05,0.12247449};

    //параметры окружающей плазмы
    dipole_physical_environment_context_t physical_environment_cntx = {
          .R0 = R0, .V0 = V0
        , .omega_cc0 = 1.0, .omega_pc0 = 0.1
        , .phi_width = 0.23, .L_width = 0.005
        , .source_density_coeff = 0.2 
    };

    /*homogeneous_physical_environment_context_t physical_environment_cntx = {
          .R0 = R0, .V0 = V0
        , .omega_cc0 = 1., .omega_pc0 = 0.1
        , .cold_density = 0.01, .source_density = 0.99
    };*/

    //контекс для дисперсионного уравнения
    /*epsilon_context_t dispersion_relation_cntx = {
        dipole_physical_environment, &physical_environment_cntx
    };*/
    epsilon_context_t dispersion_relation_cntx = {
        homogeneous_physical_environment, &physical_environment_cntx
    };


    double * points = NULL; unsigned count = 0;

    omega_correctorH_context_t corrector_cntx = {
          1.e-9, 1.e-9, 1.e-12, 1000u
        , warm_dispersion_relationH, &dispersion_relation_cntx
    };

    VectorSp R = {1.959273, 0.447126, 0.004595};

    points = calloc(4000000u, 2*sizeof(double));

    FILE * fd = fopen("./map.txt","w");
    VectorH K = {1.,1.}; double kmod = hypot(K.pl,K.pr);
    K.pl /= kmod; K.pr /= kmod;

    for (double k = 1.e-3; k < 2.; k += 2.5e-5) {
        for (double w = 0.; w < 1.5; w += 5.e-3) {
            double w_vec[2] = {w,0.};
            VectorH Kh = {K.pl*k,K.pr*k};
            if (omega_correctorH(&corrector_cntx,&R,&Kh,w_vec)) {
                if (w_vec[0] > 0 && !linear_search(k,w_vec[0],points,count,1.e-4)) {
                    fprintf(fd,"%.8f %.8f %.8f\n",k,w_vec[0],hypot(Kh.pl,Kh.pr)/w_vec[0]);
                    points[2*count] = k; points[2*count+1] = w_vec[0];
                    ++count;
                }
            }
        }
    }
    fclose(fd);
    free(points);



    //kpl_start <--> kpl_stop 
    /*****************  //w_stop
     *               *
     *               *
     *               *
     *               *  //w_start
     * ***************/
    /*struct map_domain domain = {-1.0,-0.,0.5,1.5,5.e-2,5.e-4};
    kpl_correctorH_context_t corrector_cntx = {
          1.e-9, 1.e-12, 1000u
        , warm_dispersion_relationH, &dispersion_relation_cntx
    };
    unsigned max_points = (unsigned)((domain.kpl_stop - domain.kpl_start)*(domain.w_stop-domain.w_start)/domain.kpl_step/domain.w_step)+2u;
    printf("%u\n",max_points); fflush(stdout);

    points = calloc(max_points, 2*sizeof(double));

    FILE * fd = fopen("./map.txt","w");
    VectorSp R = {2.07341,0.4491,0.11383};//R0;
    VectorH K = {-0.73224,hypot(0.61481,-0.1311)};//{0.,hypot(0.4,0.)};
    bool proceed = true;
    for (double kpl = domain.kpl_start; kpl < domain.kpl_stop && proceed; kpl += domain.kpl_step) {
        for (double w = domain.w_start; w < domain.w_stop && proceed; w += domain.w_step) {
            double w_vec[2] = {w,0.};
            VectorH Kh = {kpl,K.pr};
            if (kpl_correctorH(&corrector_cntx,&R,&Kh,w_vec)) {
                if (!linear_search(kpl,w,points,count,1.e-4)) {
                    fprintf(fd,"%.8f %.8f\n",Kh.pl,w);
                    points[2*count] = Kh.pl; points[2*count+1] = w;
                    ++count;
                }
            }
            if (count == max_points) {
                proceed = false;
                printf("Not enough memory\n");
            }
        }
    }
    fclose(fd);
    free(points);*/

    printf("%u\n",count);

    return 0;
}