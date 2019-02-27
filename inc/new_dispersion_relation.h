#ifndef IKI_DispersionRelation_H
#define IKI_DispersionRelation_H

/**
 * Asking for a real frequency (omega)
 * And returning a squared refractive index (Npr*Npr)
*/
double warm_plus_Npr(PhysicalEnvironment env, VectorH K, double w);
double warm_minus_Npr(PhysicalEnvironment env, VectorH K, double w);

#endif


