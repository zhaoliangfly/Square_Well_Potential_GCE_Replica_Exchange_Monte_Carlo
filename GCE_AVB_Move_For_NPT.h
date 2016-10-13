#ifndef GCE_AVB_Move_For_NPT_H_
#define GCE_AVB_Move_For_NPT_H_

#include "stdio.h"
#include "stdlib.h"
#include "time.h"
#include "math.h"
#include "string.h"
#include "limits.h"

#include "Parameters.h"
#include "Distance.h"
#include "dSFMT.h"
#include "Output_Energy_Configuration.h"

void GCE_AVB_Move_For_NPT_Function (int *steps, Particle_Model *p_head, dsfmt_t *dsfmt, double *dr, double *accept, double *reject, double *bx, double *by, double *bz)
{
	if (*steps > INT_MAX)
    {
        printf ("ERROR! Step exceeds INT_MAX!\n");
        exit (EXIT_SUCCESS);
    }

    int selected, i, k;
    selected = (int)(dsfmt_genrand_open_open (dsfmt) * Particle_Number);

    int i, j;
    Particle_Model temp[Particle_Number];
    for (i = 0; i < Particle_Number; i++)
    	for (j = 0; j < 3; j++)
    		temp[i].Position[j] = p_head[i].Position[j] * (boxn / (*bx));

    /**Judge if the Spheres are too close! **/
	int l = 0;
    for ( i = 0; i < Particle_Number; i++)
    {
        for ( j = i + 1; j < Particle_Number; j++)
        {
            	if (Distance_3D_1 (temp[i].Position[0], temp[i].Position[1], temp[i].Position[2], temp[j].Position[0], temp[j].Position[1], temp[j].Position[2], boxn, boxn, boxn) < 2.0 * Sigma_Sphere)
            	{
            		l++;
            		if (l != 0)
            			break;
            	}	
        }
    }

    if (l != 0) /** Two spheres are too close ! **/
    {
    	former_energy += 0.0;
        double temp_enthalpy;
        temp_enthalpy = former_energy + P_Reference * v_old;  
        Output_Energy_To_File_or_Screen (0, Output_Interval_Steps, *steps, temp_enthalpy);
        //if ((*steps > Statistic_start_step) && (*steps < Statistic_end_step))
        (*reject) += 1.0;
        (*steps) += 1;
    }
    else
    {
    	double rand_g, acc;
    	rand_g = dsfmt_genrand_open_open (dsfmt);
    	double delta_energy, e_old, e_new, delta_enthalpy;
    	e_old = former_energy;
    	e_new = Total_Energy_Square_Well_Potential (&(temp[0]), boxn, boxn, boxn);
    	delta_energy = e_new - e_old;
    	delta_enthalpy = delta_energy + P_Reference * (v_new - v_old);
    	/** ATTENTION!!!!!!!!---------- GCE!!!!!! **/
        //acc = exp (-1.0 * delta_energy * BETA_RT - 0.5 * Alpha * delta_energy * (delta_energy + 2.0 * former_energy - 2.0 * E0));
        /** End **/

        /** ATTEMTION!!!!!!!!---------- GCE in NPT ensemble! **/
        acc = exp (-1.0 * delta_enthalpy * BETA_RT - 0.5 * Alpha * delta_enthalpy * (delta_enthalpy + 2.0 * (former_energy + P_Reference * v_old) - 2.0 * H0) - (Particle_Number + 1) * log (v_new / v_old) / BETA_RT);
        /** End **/
    	if (rand_g < acc) /** Acceptance **/
    	{
    		former_energy = e_new;
            double temp_enthalpy_1;
            temp_enthalpy_1 = former_energy + P_Reference * v_new; 
    		Output_Energy_To_File_or_Screen (0, Output_Interval_Steps, *steps, temp_enthalpy_1);
    		//if ((*steps > Statistic_start_step) && (*steps < Statistic_end_step))
    		(*accept) += 1.0;
    		(*steps) += 1;

    		for (i = 0; i < Particle_Number; i++)
    			for (j = 0; j < 3; j++)
    				p_head[i].Position[j] = temp[i].Position[j];
    		
    		*bx = boxn;
    		*by = boxn;
    		*bz = boxn;
        }
    	else
    	{
    		former_energy += 0.0;
            double temp_enthalpy_2;
            temp_enthalpy_2 = former_energy + P_Reference * v_old;
    		Output_Energy_To_File_or_Screen (0, Output_Interval_Steps, *steps, temp_enthalpy_2);
    		//if ((*steps > Statistic_start_step) && (*steps < Statistic_end_step))
    		(*reject) += 1.0;
    		(*steps) += 1;	
    	}
    }

}




#endif