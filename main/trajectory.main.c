#include <physical_environment.h>
#include <dispersion_relation.h>
#include <correction.h>
#include <ray_trace.h>
#include <vector.h>

#include <stdio.h>
#include <math.h>

static inline
double sqr(double x) { return x*x; }

VectorSp Ksp_from_KqL(VectorSp H, VectorH KqL, double Kphi) {
    double Hmod = hypot(H.r,H.th);
    VectorSp K_Sp_pl = {H.r*KqL.pl/Hmod, H.th*KqL.pl/Hmod, 0.};
    VectorSp K_Sp_pr = {H.th*KqL.pr/Hmod, -H.r*KqL.pr/Hmod, 0.};
    return (VectorSp){K_Sp_pl.r+K_Sp_pr.r, K_Sp_pl.th+K_Sp_pr.th, Kphi};
}

int main(int argc, char **args) {
    //ФИЗИЧЕСКАЯ МОДЕЛЬ
    VectorSp R0 = {1.975,0.44928364,0.}; //точка нормировки
    VectorH V0  = {0.05,0.12247449};     //скорость пучка в точке нормировки

    //параметры окружающей плазмы и каверны
    dipole_physical_environment_context_t physical_environment_cntx = {
          .R0 = R0, .V0 = V0
        , .omega_cc0 = 1.0, .omega_pc0 = 0.1
        , .phi_width = 0.23, .L_width = 0.005
        , .source_density_coeff = 0.2 
    };

    //контекс для дисперсионного уравнения
    epsilon_context_t dispersion_relation_cntx = {
        dipole_physical_environment, &physical_environment_cntx
    };


    double r, th, phi, kq, kL, kphi, w, g; //переменные для получения данных
    //получаем значения точки, с которой необходимо начать вычисления
    {
        FILE *initial_points_stream = NULL;
        if (argc > 1) {
            initial_points_stream = fopen(args[1],"r");
        }

        if (NULL == initial_points_stream) {
            fprintf(stdout,"Can't open file!\n Reading data from standard input\n");
            initial_points_stream = stdout;
        }
        
        fscanf(initial_points_stream,"%lf %lf %lf %lf %lf %lf %lf %lf",&r,&th,&phi,&kq,&kL,&kphi,&w,&g);
    
        if (initial_points_stream != stdout) {
            fclose(initial_points_stream);
        }
    }
    VectorSp R = {r,th,phi};
    VectorH Kh = {kq,hypot(kL,kphi)};
    double w_vec[2] = {w,g};

    PhysicalEnvironment starting_point_environment = dipole_physical_environment(&physical_environment_cntx,R);
    VectorSp K = Ksp_from_KqL(starting_point_environment.H,(VectorH){kq,kL},kphi);

    //проверка введённого значения комплексной частоты
    {
        omega_correctorH_context_t omega_correctorH_cntx = {
            1.e-9,1.e-9,1.e-12,1000u
            , warm_dispersion_relationH, &dispersion_relation_cntx
        };
        
        if (!omega_correctorH(&omega_correctorH_cntx,&R,&Kh,w_vec)) {
            printf("Terminated! Bad initial guess: %.9f %.9f %.9f => %.9f %.9f\n",kq,kL,kphi,w,g);
            return 0;
        }

        printf(
              "(r,th,phi) = {%.9f, %.9f, %.9f}\n(kq,kL,kphi) = {%.9f,%.9f,%.9f}\n(kr,kth,kphi) = {%.9f,%.9f,%.9f}\nw = %.9f + %.9f*i\n"
            , r,th,phi, kq,kL,kphi, K.r,K.th,K.phi, w_vec[0], w_vec[1]
        );
        fflush(stdout);
        
        char buff[2];
        printf("Do you want to proceed? "); fflush(stdout);
        scanf("%1s",buff);
        if ('Y' != buff[0] && 'y' != buff[0]) return 0;
    }
    
    double t = 0., S = 0.; //время и коэффициент усиления общие между траекториями

    //Рассчёт участка траектории с набором энергии
    FILE *gain_trj = fopen("./gain.txt","w");
    if (NULL == gain_trj) {
        printf("Can't open file to write gain data\n");
        return 0;
    }

    fprintf(gain_trj,"t dt r th phi q L Kr Kth Kphi Kq KL gamma S Vr Vth Vphi VKr VKth VKphi d_cold d_source N\n");
    {
        VectorSp VR, VK; double dt = 2.5e-5, g = w_vec[1];

        //корректор k-параллельного и инкремента контекст
        gamma_kpr_correctorSp_context_t gain_corrector_cntx = {
            dipole_physical_environment, &physical_environment_cntx
            ,  1.e-9, 1.e-9, 1.e-12, 1000u
            , warm_dispersion_relationH, &dispersion_relation_cntx
        };
        
        //предиктор-корректор контекст
        predictor_step_context_t predictor_step_cntx = {
              (VectorSp){1.e-6,1.e-6,1.e-6}, (VectorSp){1.e-6,1.e-6,1.e-6}, 1.e-6
            , warm_dispersion_relationSp, &dispersion_relation_cntx
            , gamma_kpr_correctorSp, &gain_corrector_cntx
        };

        
        while(w_vec[1] > 1.e-4) {
            
            step_status_t status = predictor_step(&predictor_step_cntx, &R,  &K, w_vec, &VR, &VK, dt);
            
            t += dt;  printf("%f %d %f\n", t, status, w_vec[1]);
            if (STP_OK != status) { printf("ERROR"); return 0; }

            S += 0.5*(w_vec[1]+g)*dt*2.*6380./3.*2.*3.14159265; g = w_vec[1];

            PhysicalEnvironment phys_env = dipole_physical_environment(&physical_environment_cntx,R);
            VectorSp H = phys_env.H;
            
            double hmod = hypot(H.r,H.th);
            VectorSp He = {H.r/hmod,H.th/hmod}, Hp = {H.th/hmod,-H.r/hmod};

            fprintf(gain_trj,"%f %f ", t, dt);
            fprintf(gain_trj,"%f %f %f %f %f ", R.r, R.th, R.phi, 0.5*sqr(1./R.r)*cos(R.th), R.r/sqr(sin(R.th)));
            fprintf(gain_trj,"%f %f %f %f %f ", K.r, K.th, K.phi, projection_of_on(K,He).pl, projection_of_on(K,Hp).pl);
            fprintf(gain_trj,"%f %f ", w_vec[1], S);
            fprintf(gain_trj,"%f %f %f %f %f %f ", VR.r ,VR.th, VR.phi, VK.r, VK.th, VK.phi);
            fprintf(gain_trj,"%f %f %f ", phys_env.density_c, phys_env.density_s, sqrt(vecsp_norm(K))/w_vec[0]);
            fprintf(gain_trj,"\n");

            if (sqrt(vecsp_norm(K))/w_vec[0] >= 1. ) {
                printf("N > 1\n"); fflush(stdout);
                scanf("%*s");
                break;
            }
        }
    }
    fclose(gain_trj);


    //Рассчёт пути без набора энергии
    FILE *stable_trj = fopen("./stable.txt","w");
    if (NULL == stable_trj) {
        printf("Can't open file to save stable trajectory\n");
        return 0;
    }

    fprintf(stable_trj,"t dt r th phi q L Kr Kth Kphi Kq KL gamma S Vr Vth Vphi VKr VKth VKphi d_cold d_source N\n");
    {
        VectorSp VR, VK; double dt = 0.5e-5; w_vec[1] = 0.;

        N_context_t n_dispersion_relation_cntx = {
              dipole_physical_environment, &physical_environment_cntx
            , sqrt(vecsp_norm(K))/w_vec[0]
        };

        kpl_correctorSp_context_t n_corrector_cntx = {
            dipole_physical_environment, &physical_environment_cntx
            ,  1.e-9, 1.e-12, 1000u
            , source_dispersion_relationH, &dispersion_relation_cntx
        };

        predictor_step_context_t n_predictor_step_cntx = {
              (VectorSp){1.e-6,1.e-6,1.e-6}, (VectorSp){1.e-6,1.e-6,1.e-6}, 1.e-6
            , source_dispersion_relationSp, &dispersion_relation_cntx
            , kpl_correctorSp, &n_corrector_cntx
        };
        
        kpl_correctorSp_context_t cold_corrector_cntx = {
              dipole_physical_environment, &physical_environment_cntx
            , 1.e-9, 1.e-12, 1000u
            , cold_dispersion_relationH, &dispersion_relation_cntx
        };

        predictor_step_context_t cold_predictor_step_cntx = {
              (VectorSp){1.e-6,1.e-6,1.e-6}, (VectorSp){1.e-6,1.e-6,1.e-6}, 1.e-6
            , cold_dispersion_relationSp, &dispersion_relation_cntx
            , kpl_correctorSp, &cold_corrector_cntx     
        };
        
        bool proceed = true;
        while (t < 1. && proceed)
        {
            step_status_t status = predictor_step(&n_predictor_step_cntx, &R,  &K, w_vec, &VR, &VK, dt);
            
            t += dt;  printf("%f %d %f\n", t, status, w_vec[1]);
            if (STP_OK != status) { printf("ERROR"); return 0; }

            PhysicalEnvironment phys_env = dipole_physical_environment(&physical_environment_cntx,R);
            VectorSp H = phys_env.H;
            
            double hmod = hypot(H.r,H.th);
            VectorSp He = {H.r/hmod,H.th/hmod}, Hp = {H.th/hmod,-H.r/hmod};

            fprintf(stable_trj,"%f %f ", t, dt);
            fprintf(stable_trj,"%f %f %f %f %f ", R.r, R.th, R.phi, 0.5*sqr(1./R.r)*cos(R.th), R.r/sqr(sin(R.th)));
            fprintf(stable_trj,"%f %f %f %f %f ", K.r, K.th, K.phi, projection_of_on(K,He).pl, projection_of_on(K,Hp).pl);
            fprintf(stable_trj,"%f %f ", w_vec[1], S);
            fprintf(stable_trj,"%f %f %f %f %f %f ", VR.r ,VR.th, VR.phi, VK.r, VK.th, VK.phi);
            fprintf(stable_trj,"%f %f %f ", phys_env.density_c, phys_env.density_s, sqrt(vecsp_norm(K))/w_vec[0]);
            fprintf(stable_trj,"\n");

            {
                VectorSp K_cold = K;
                if (kpl_correctorSp(&cold_corrector_cntx,&R,&K_cold,w_vec) && fabs(K.r-K_cold.r) < 5.e-4) {
                    printf("K     = %f %f %f\n",K.r,K.th,K.phi);
                    printf("Kcold = %f %f %f\n",K_cold.r,K_cold.th,K_cold.phi);
                    fflush(stdout);
                    {
                        char buffer[2];
                        scanf("%1s",buffer);
                        if ('y' == buffer[0] || 'Y' == buffer[0]) {
                            fprintf(stable_trj,"\n");
                            K = K_cold;
                            proceed = false;
                        }
                    }
                }
            }
        }

        //dt = 2.5e-5;
        while (t < 1.)
        {
            step_status_t status = predictor_step(&cold_predictor_step_cntx, &R,  &K, w_vec, &VR, &VK, dt);
            
            t += dt;  printf("%f %d %f\n", t, status, w_vec[1]);
            if (STP_OK != status) { printf("ERROR"); return 0; }

            PhysicalEnvironment phys_env = dipole_physical_environment(&physical_environment_cntx,R);
            VectorSp H = phys_env.H;
            
            double hmod = hypot(H.r,H.th);
            VectorSp He = {H.r/hmod,H.th/hmod}, Hp = {H.th/hmod,-H.r/hmod};

            fprintf(stable_trj,"%f %f ", t, dt);
            fprintf(stable_trj,"%f %f %f %f %f ", R.r, R.th, R.phi, 0.5*sqr(1./R.r)*cos(R.th), R.r/sqr(sin(R.th)));
            fprintf(stable_trj,"%f %f %f %f %f ", K.r, K.th, K.phi, projection_of_on(K,He).pl, projection_of_on(K,Hp).pl);
            fprintf(stable_trj,"%f %f ", w_vec[1], S);
            fprintf(stable_trj,"%f %f %f %f %f %f ", VR.r ,VR.th, VR.phi, VK.r, VK.th, VK.phi);
            fprintf(stable_trj,"%f %f %f ", phys_env.density_c, phys_env.density_s, sqrt(vecsp_norm(K))/w_vec[0]);
            fprintf(stable_trj,"\n");
        }
    }

    fclose(stable_trj);

    return 0;
}