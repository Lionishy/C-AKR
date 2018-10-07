#ifndef IKI_PhysicalEnvironment_H
#define IKI_PhysicalEnvironment_H

#include <vector.h>

typedef struct _PhysicalEnvironment {
    VectorSp H;
    double cavity, density_c, density_s;
    double omega_cc, omega_pc, omega_ps;
    VectorH V;
} PhysicalEnvironment;

typedef PhysicalEnvironment (*phys_env_t)(void const *cntx, VectorSp R);


typedef struct _dipole_physical_environment_context {
    VectorSp R0;
    VectorH V0;
    double omega_cc0, omega_pc0;
    double phi_width, L_width, source_density_coeff;
} dipole_physical_environment_context_t;

PhysicalEnvironment dipole_physical_environment(void const *cntx, VectorSp R);


typedef struct _homogeneous_physical_environment_context {
    VectorSp R0;
    VectorH V0;
    double omega_cc0, omega_pc0;
    double cold_density, source_density;
} homogeneous_physical_environment_context_t;

PhysicalEnvironment homogeneous_physical_environment(void const *cntx, VectorSp R);

#endif // IKI_PhysicalEnvironment_H