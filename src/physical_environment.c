#include <physical_environment.h>
#include <vector.h>

#include <math.h>

PhysicalEnvironment dipole_physical_environment(void const * ptr, VectorSp R) {
    dipole_physical_environment_context_t const * cntx = ptr;
    
    VectorSp H = {-cos(R.th)/(R.r*R.r*R.r), -0.5*sin(R.th)/(R.r*R.r*R.r), 0.};
    VectorSp H0 = {-cos(cntx->R0.th)/(cntx->R0.r*cntx->R0.r*cntx->R0.r), -0.5*sin(cntx->R0.th)/(cntx->R0.r*cntx->R0.r*cntx->R0.r), 0.};
    
    double L = acos(sin(R.th)*sqrt(1./R.r)), L0 = acos(sin(cntx->R0.th)*sqrt(1./cntx->R0.r));
    double cavity = exp( -(L-L0)*(L-L0)/(cntx->L_width*cntx->L_width) -(R.phi-cntx->R0.phi)*(R.phi-cntx->R0.phi)/(cntx->phi_width*cntx->phi_width) );
    double 
        density_c = cntx->R0.r*cntx->R0.r/(R.r*R.r)*(1.-cavity),
        density_s = cntx->R0.r*cntx->R0.r/(R.r*R.r)*cntx->source_density_coeff*cavity;

    double 
        omega_cc = cntx->omega_cc0*hypot(H.r,H.th)/hypot(H0.r,H0.th),
        omega_pc = cntx->omega_pc0*sqrt(density_c),
        omega_ps = cntx->omega_pc0*sqrt(density_s);

    VectorH V;
    //velocity squared !!!
    V.pr = cntx->V0.pr*cntx->V0.pr*hypot(H.r,H.th)/hypot(H0.r,H0.th);
    V.pl = cntx->V0.pr*cntx->V0.pr + cntx->V0.pl*cntx->V0.pl - V.pr;
    V.pr = sqrt(V.pr); V.pl = sqrt(V.pl);

    return (PhysicalEnvironment){H,cavity,density_c,density_s,omega_cc,omega_pc,omega_ps,V};
}

PhysicalEnvironment homogeneous_physical_environment(void const *ptr, VectorSp R) {
    homogeneous_physical_environment_context_t const * cntx = ptr;
    VectorSp H0 = {-cos(cntx->R0.th)/(cntx->R0.r*cntx->R0.r*cntx->R0.r), -0.5*sin(cntx->R0.th)/(cntx->R0.r*cntx->R0.r*cntx->R0.r), 0.};
    double 
        omega_pc = cntx->omega_pc0*sqrt(cntx->cold_density),
        omega_ps = cntx->omega_pc0*sqrt(cntx->source_density);
    return (PhysicalEnvironment){H0,1.,cntx->cold_density,cntx->source_density,cntx->omega_cc0,omega_pc,omega_ps,cntx->V0};
}