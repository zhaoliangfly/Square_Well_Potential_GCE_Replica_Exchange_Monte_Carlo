#include "math.h"
#include "stdio.h"

#include "Parameters.h"
#include "Distance.h"

#ifndef Square_Well_Potential_Definition_H_
#define Square_Well_Potential_Definition_H_

/** The first realization of GCE is based on the Square well potential model! **/

double Square_Well_Potential(Particle_Model *p_i, Particle_Model *p_j, double bx, double by, double bz)
{
	double d;
	d = Distance_3D_1 (p_i -> Position[0], p_i -> Position[1], p_i -> Position[2], p_j -> Position[0], p_j -> Position[1], p_j -> Position[2], bx, by, bz);
	if ((d > 2.0 * Sigma_Sphere) && (d < F_Sigma_Sphere * 2.0 * Sigma_Sphere))
		return -1.0 * Epsilon;
	if ( d > F_Sigma_Sphere * 2.0 * Sigma_Sphere)
		return 0.0; 
}

double Delta_Energy_Single_Particle_Square_Well_Potential(int selected_number, Particle_Model *p_head, Particle_Model *p_new, double bx, double by, double bz)
{
	double e_old, e_new, delta_e;
	e_old = 0.0;
	e_new = 0.0;
	int i;
	for (i = 0; i < Particle_Number; i++)
	{
		if (i != selected_number)
		{
			e_old += Square_Well_Potential(p_head + selected_number, p_head + i, bx, by, bz);
			e_new += Square_Well_Potential(p_new, p_head + i, bx, by, bz); 
		}
	}
	delta_e = e_new - e_old;
	return delta_e;
}

double Total_Energy_Square_Well_Potential (Particle_Model *p_head, double bx, double by, double bz)
{
	int i, j;
	double e = 0.0;
	for (i = 0; i < Particle_Number; i++)
		for (j = i + 1; j < Particle_Number; j++)
			e += Square_Well_Potential(p_head + i, p_head + j, bx, by, bz);
	return e;
}

#endif
