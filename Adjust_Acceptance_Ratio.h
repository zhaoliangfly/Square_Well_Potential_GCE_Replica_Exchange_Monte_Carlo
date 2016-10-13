#ifndef Adjust_Acceptance_Ratio_H_
#define Adjust_Acceptance_Ratio_H_

#include "stdio.h"
#include "math.h"
#include "Parameters.h"

void Adjust_Acceptance (int steps, int type, double *accept, double *reject, double *dr)
{
        if (steps % Adjust_Ratio_Interval_Steps == 0)
        {
                double current_ratio;
                current_ratio = (*accept)/(*accept + *reject);
                if (current_ratio > Object_Ratio)
                        (*dr) *= 1.05;  //length of the translation, unit nm.
                if (current_ratio < Object_Ratio)
                        (*dr) *= 0.95;
                if (*dr < 0.01) // The minimum is 0.02 nm.
                {
                        //printf("WARNING!");
                        (*dr) = 0.01;
                }
                if (*dr > 2.5)
                {
                        (*dr) = 2.5;
                }

                if (type == 0)
                {
                        FILE *fp;
                        char str_head[80] = "Result_Acceptance_Ratio";
                        char al[15];
                        sprintf (al, "%.3f", Alpha);
                        char be[15];
                        sprintf (be, "%.3f", BETA_RT);
                        char pressure[15];
                        sprintf (pressure, "%.3f", P_Reference);
                        char e[15];
                        sprintf (e, "%.3f", fabs (H0));
                        strcat (str_head, "_AL_");
                        strcat (str_head, al);
                        strcat (str_head, "_BE_");
                        strcat (str_head, be);
                        strcat (str_head, "_Pr_");
                        strcat (str_head, pressure);
                        strcat (str_head, "_H0_");
                        strcat (str_head, e);
                        strcat (str_head, ".dat");

                        fp = fopen(str_head, "a");
                        fprintf(fp, "%d %lf %lf\n", steps, current_ratio, *dr);
                        fclose (fp);
                }

                if (type == 1)
                {
                        printf("%d %lf %lf\n", steps, current_ratio, *dr);
                }

                if (type == 2)
                {
                        //No operation.
                }

                *accept = 0;
                *reject = 0;
        }
}


void Adjust_Acceptance_Mpi (int id, int type, int interval_steps, int steps, double *accept, double *reject, double *dr)
{
        if (steps % interval_steps == 0)
        {
                double current_ratio;
                current_ratio = (*accept)/(*accept + *reject);
                if (current_ratio > Object_Ratio)
                        (*dr) *= 1.05;  //length of the translation, unit nm.
                if (current_ratio < Object_Ratio)
                        (*dr) *= 0.95;
                if (*dr < 0.01) // The minimum is 0.02 nm.
                {
                        //printf("WARNING!");
                        (*dr) = 0.01;
                }
                if (*dr > 2.5)
                {
                        (*dr) = 2.5;
                }

                if (type == 0)
                {
                        FILE *fp;
                        char str_head[80] = "Result_Acceptance_Ratio";
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

                        fp = fopen(str_head, "a");
                        fprintf(fp, "%d %lf %lf\n", steps, current_ratio, *dr);
                        fclose (fp);
                }

                if (type == 1)
                {
                        printf("%d %lf %lf\n", steps, current_ratio, *dr);
                }

                if (type == 2)
                {
                        //No operation.
                }

                *accept = 0;
                *reject = 0;
        }
}


#endif
