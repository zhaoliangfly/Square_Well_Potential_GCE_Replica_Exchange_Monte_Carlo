#ifndef GCE_Metropolis_Move_For_NPT_H_
#define GCE_Metropolis_Move_For_NPT_H_

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

void GCE_Metropolis_Move_For_NPT_Function (int *steps, Particle_Model *p_head, dsfmt_t *dsfmt, double *dr, double *accept, double *reject, double bx, double by, double bz)
{
    if (*steps > INT_MAX)
    {
        printf ("ERROR! Step exceeds INT_MAX!\n");
        exit (EXIT_SUCCESS);
    }

	int selected, i, k;
	selected = (int)(dsfmt_genrand_open_open (dsfmt) * Particle_Number);

	Particle_Model temp;
	for (i = 0; i < 3; i++)
		temp.Position[i] = (p_head + selected) -> Position[i] + (dsfmt_genrand_open_open (dsfmt) - 0.5) * (*dr);

	/* Check the PBC */
	if (temp.Position[0] > bx) temp.Position[0] -= bx; /**Check the Periodical Boundary Condition**/
	if (temp.Position[0] < 0.0) temp.Position[0] += bx;
	if (temp.Position[1] > by) temp.Position[1] -= by;
	if (temp.Position[1] < 0.0) temp.Position[1] += by;
	if (temp.Position[2] > bz) temp.Position[2] -= bz;
	if (temp.Position[2] < 0.0) temp.Position[2] += bz;

	/**Judge if the Spheres are too close! **/
	int l = 0;
    for ( i = 0; i < Particle_Number; i++)
    {
        if ( i != selected)
        {
            	if (Distance_3D_1 (temp.Position[0], temp.Position[1], temp.Position[2], (p_head + i) -> Position[0], (p_head + i) -> Position[1], (p_head + i) -> Position[2], bx, by, bz) < 2.0 * Sigma_Sphere)
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
        temp_enthalpy = former_energy + P_Reference  * bx * by * bz;  
        Output_Energy_To_File_or_Screen (To_File_Or_Screen, Output_Interval_Steps, *steps, temp_enthalpy);
        //if ((*steps > Statistic_start_step) && (*steps < Statistic_end_step))
        (*reject) += 1.0;
        (*steps) += 1;
    }
    else
    {
    	double rand_g, acc;
    	rand_g = dsfmt_genrand_open_open (dsfmt);
    	double delta_energy;
    	delta_energy = Delta_Energy_Single_Particle_Square_Well_Potential (selected, p_head, &temp, bx, by, bz);
    	/** ATTENTION!!!!!!!!---------- GCE!!!!!! **/
        //acc = exp (-1.0 * delta_energy * BETA_RT - 0.5 * Alpha * delta_energy * (delta_energy + 2.0 * former_energy - 2.0 * E0));
        /** End **/

        /** ATTEMTION!!!!!!!!---------- GCE in NPT ensemble! **/
        acc = exp (-1.0 * delta_energy * BETA_RT - 0.5 * Alpha * delta_energy * (delta_energy + 2.0 * (former_energy + P_Reference * bx * by * bz) - 2.0 * H0));
        /** End **/
    	if (rand_g < acc) /** Acceptance **/
    	{
    		former_energy += delta_energy;
            double temp_enthalpy_1;
            temp_enthalpy_1 = former_energy + P_Reference * bx * by * bz; 
    		Output_Energy_To_File_or_Screen (To_File_Or_Screen, Output_Interval_Steps, *steps, temp_enthalpy_1);
    		//if ((*steps > Statistic_start_step) && (*steps < Statistic_end_step))
    		(*accept) += 1.0;
    		(*steps) += 1;

    		int m;
    		for (m = 0; m < 3; m++)
    			(p_head + selected) -> Position[m] = temp.Position[m];
    	}
    	else
    	{
    		former_energy += 0.0;
            double temp_enthalpy_2;
            temp_enthalpy_2 = former_energy + P_Reference * bx * by * bz;
    		Output_Energy_To_File_or_Screen (To_File_Or_Screen, Output_Interval_Steps, *steps, temp_enthalpy_2);
    		//if ((*steps > Statistic_start_step) && (*steps < Statistic_end_step))
    		(*reject) += 1.0;
    		(*steps) += 1;	
    	}
    }
}

void GCE_Metropolis_Move_For_NPT_Function_Mpi (int id, int *steps, Particle_Model *p_head, dsfmt_t *dsfmt, double *dr, double *accept, double *reject, double bx, double by, double bz)
{
    if (*steps > INT_MAX)
    {
        printf ("ERROR! Step exceeds INT_MAX!\n");
        exit (EXIT_SUCCESS);
    }

    int selected, i, k;
    selected = (int)(dsfmt_genrand_open_open (dsfmt) * Particle_Number);

    Particle_Model temp;
    for (i = 0; i < 3; i++)
        temp.Position[i] = (p_head + selected) -> Position[i] + (dsfmt_genrand_open_open (dsfmt) - 0.5) * (*dr);

    /* Check the PBC */
    if (temp.Position[0] > bx) temp.Position[0] -= bx; /**Check the Periodical Boundary Condition**/
    if (temp.Position[0] < 0.0) temp.Position[0] += bx;
    if (temp.Position[1] > by) temp.Position[1] -= by;
    if (temp.Position[1] < 0.0) temp.Position[1] += by;
    if (temp.Position[2] > bz) temp.Position[2] -= bz;
    if (temp.Position[2] < 0.0) temp.Position[2] += bz;

    /**Judge if the Spheres are too close! **/
    int l = 0;
    for ( i = 0; i < Particle_Number; i++)
    {
        if ( i != selected)
        {
                if (Distance_3D_1 (temp.Position[0], temp.Position[1], temp.Position[2], (p_head + i) -> Position[0], (p_head + i) -> Position[1], (p_head + i) -> Position[2], bx, by, bz) < 2.0 * Sigma_Sphere)
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
        temp_enthalpy = former_energy + P_Reference  * bx * by * bz;  
        //Output_Energy_To_File_or_Screen (2, Output_Interval_Steps, *steps, temp_enthalpy);
        Output_Energy_To_File_or_Screen_Mpi (id, To_File_Or_Screen_Mpi_Enthalpy, Output_Interval_Steps, *steps, temp_enthalpy);
        //if ((*steps > Statistic_start_step) && (*steps < Statistic_end_step))
        (*reject) += 1.0;
        (*steps) += 1;
    }
    else
    {
        double rand_g, acc;
        rand_g = dsfmt_genrand_open_open (dsfmt);
        double delta_energy;
        delta_energy = Delta_Energy_Single_Particle_Square_Well_Potential (selected, p_head, &temp, bx, by, bz);
        /** ATTENTION!!!!!!!!---------- GCE!!!!!! **/
        //acc = exp (-1.0 * delta_energy * BETA_RT - 0.5 * Alpha * delta_energy * (delta_energy + 2.0 * former_energy - 2.0 * E0));
        /** End **/

        /** ATTEMTION!!!!!!!!---------- GCE in NPT ensemble! **/
        acc = exp (-1.0 * delta_energy * BETA_RT - 0.5 * Alpha * delta_energy * (delta_energy + 2.0 * (former_energy + P_Reference * bx * by * bz) - 2.0 * H0));
        /** End **/
        if (rand_g < acc) /** Acceptance **/
        {
            former_energy += delta_energy;
            double temp_enthalpy_1;
            temp_enthalpy_1 = former_energy + P_Reference * bx * by * bz; 
            //Output_Energy_To_File_or_Screen (2, Output_Interval_Steps, *steps, temp_enthalpy_1);
            Output_Energy_To_File_or_Screen_Mpi (id, To_File_Or_Screen_Mpi_Enthalpy, Output_Interval_Steps, *steps, temp_enthalpy_1);
            //if ((*steps > Statistic_start_step) && (*steps < Statistic_end_step))
            (*accept) += 1.0;
            (*steps) += 1;

            int m;
            for (m = 0; m < 3; m++)
                (p_head + selected) -> Position[m] = temp.Position[m];
        }
        else
        {
            former_energy += 0.0;
            double temp_enthalpy_2;
            temp_enthalpy_2 = former_energy + P_Reference * bx * by * bz;
            //Output_Energy_To_File_or_Screen (2, Output_Interval_Steps, *steps, temp_enthalpy_2);
            Output_Energy_To_File_or_Screen_Mpi (id, To_File_Or_Screen_Mpi_Enthalpy, Output_Interval_Steps, *steps, temp_enthalpy_2);
            //if ((*steps > Statistic_start_step) && (*steps < Statistic_end_step))
            (*reject) += 1.0;
            (*steps) += 1;  
        }
    }
}


#endif
