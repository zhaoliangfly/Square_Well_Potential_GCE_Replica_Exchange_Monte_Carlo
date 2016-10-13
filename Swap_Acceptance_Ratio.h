#ifndef Swap_Acceptance_Ratio_H_
#define Swap_Acceptance_Ratio_H_

#include "stdio.h"
#include "math.h"
#include "Parameters.h"

void Output_Swap_Acceptance_Ratio_Mpi(int id, int numproc, int type, int steps, double *accept_swap, double *reject_swap)
{
	if (steps % Calculate_Swap_Acceptance_Ratio_Interval_Steps == 0)
	{
		double current_ratio_swap;
		current_ratio_swap = (*accept_swap)/(*accept_swap + *reject_swap);
		if (type == 0)
		{
		    if (id < numproc - 1)
	            {
                        FILE *fp;
                        char str_head[100] = "Result_Swap_Acceptance_Ratio";
                        char al[15];
                        sprintf (al, "%.3f", Alpha);
                        char be[15];
                        sprintf (be, "%.3f", BETA_RT);
                        char pressure[15];
                        sprintf (pressure, "%.3f", P_Reference);
                        char e[15];
                        sprintf (e, "%.3f", fabs (H0));
                        char mpiid[15];
                        sprintf (mpiid, "%d", id);
                        strcat (str_head, "_AL_");
                        strcat (str_head, al);
                        strcat (str_head, "_BE_");
                        strcat (str_head, be);
                        strcat (str_head, "_Pr_");
                        strcat (str_head, pressure);
                        strcat (str_head, "_H0_");
                        strcat (str_head, e);
                        strcat (str_head, "_Mpi_id_");
                        strcat (str_head, mpiid);
                        strcat (str_head, ".dat");

                        fp = fopen (str_head, "a");
			fprintf (fp, "%d %lf\n", steps, current_ratio_swap);
                        fclose (fp);
                    }
		}
		if (type == 1)
		{
			//printf("%d %d %lf\n", id, steps, current_ratio_swap);
			if (id < numproc - 1)
			{
				printf("%d %d %lf\n", id, steps, current_ratio_swap);
			}			
		}
		if (type == 2)
		{
			//No operation!
		}
		*accept_swap = 0;
        *reject_swap = 0;
	}
}

#endif
