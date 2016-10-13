#ifndef Initialization_System_H_
#define Initialization_System_H_

#include "stdio.h"
#include "math.h"
#include "limits.h"

#include "dSFMT.h"
#include "Parameters.h"
#include "Distance.h"
#include "Output_Energy_Configuration.h"

void Initialize_Particle_Model (int type, Particle_Model *p_head, double bx, double by, double bz)
{
	if (type == 0) /** Random Distribution **/
	{
		int i,j,k;
	    for (i = 0; i < Particle_Number; i++)
	    	for (j = 0; j < Coordinate_Number; j++)
				(p_head + i) -> Position[j] = 0.0;
	
	/* First to Srand the Position of the Particles *
	 * The x,y,z are located at the interval of     *
	 * Box_X, Box_Y, Box_Z.                         *
	 */
	 	i = 1;
		double b[3]; 
		b[0] = bx; b[1] = by; b[2] = bz;
    	dsfmt_t dsfmt; /** Predefinition to generate a random number! **/
		dsfmt_init_gen_rand (&dsfmt, (unsigned)time(NULL));
		while (i <= Particle_Number)
		{
			if (i == 1)
			{
				for (k = 0; k < 3; k++)
				{
					(p_head + i - 1) -> Position[k] = (dsfmt_genrand_open_open (&dsfmt) * b[k]);
				}
				i++;
			}
			else
			{
				int seed_time = 0;
				while (1)
				{
					if (seed_time == Initial_Position_Max)
					{
						printf ("Sorry, Cannot Fill the Box!!!\n");
						exit (EXIT_SUCCESS);
					}
					double temp[3];
					for(k = 0; k < 3; k++)
					{
						temp[k] = (dsfmt_genrand_open_open (&dsfmt) * b[k]);
			    	}
					int success = 0;
					for (j = 1;j < i; j++)
					{
						if (Distance_3D_1 (temp[0], temp[1], temp[2], (p_head + j - 1) -> Position[0], (p_head + j - 1) -> Position[1], (p_head + j - 1) -> Position[2], bx, by, bz) > 2.0 * Sigma_Sphere)
							success ++;
						else
						{
							seed_time ++;
							break;
						}
					}
					if (success == i-1)
					{
						(p_head + i - 1) -> Position[0] = temp[0];
						(p_head + i - 1) -> Position[1] = temp[1];
						(p_head + i - 1) -> Position[2] = temp[2];
						i++;
						break;
					}
				}
			}
		}

		//Output_Output_Configuration_To_File (0, 0, p_head, bx, by, bz);
	}

	if (type == 1)
	{
		FILE *fp_read;
		fp_read = fopen ("Initial_For_Read.gro","r");
		if (fp_read == NULL)
		{
			printf ("ERROR! No Initial Configuration Input!\n");
			exit (EXIT_SUCCESS);
		}
		int l;
		for (l = 0; l < Particle_Number; l++)
			fscanf (fp_read, "%lf %lf %lf\n", &((p_head + l) -> Position[0]), &((p_head + l) -> Position[1]), &((p_head + l) -> Position[2]));
		fclose (fp_read);
	}   
}

#endif
