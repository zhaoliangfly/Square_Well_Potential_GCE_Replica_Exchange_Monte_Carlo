#ifndef Output_Energy_Configuration_H_
#define Output_Energy_Configuration_H_

#include "stdio.h"
#include "math.h"
#include "stdlib.h"
#include "string.h"
#include "limits.h"

#include "Parameters.h"
#include "Square_Well_Potential_Definition.h"

void Output_Energy_To_File_or_Screen (int type, int interval_steps, int steps, double energy)
{
    /** Type 0: to File!   **
     ** Type 1: to Screen! **
     ** Type 2: Without Output **
     ** Interval_step: Output the Energy per interval steps **
     **/
    if (type == 0)
    {
        if (steps % interval_steps == 0)
        {
            char str_head[80] = "Result_Total_Enthalpy_Output";
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
            FILE *fp_energy;
            fp_energy = fopen (str_head, "a");
            fprintf (fp_energy, "%d %lf\n", steps, energy);
            fclose (fp_energy);
        }
        
    }

    if (type == 1)
    {
        if (steps % interval_steps == 0)
            printf ("%d %lf\n",steps, energy);
    }
    if (type == 2)
    {
        //Do nothing
    }
}

void Output_Radius_Distribution_Function_To_File (int steps, double *sta_rdf, Particle_Model *p_head, double bx, double by, double bz)
{
	if ((steps >= Statistic_start_step) && (steps <= Statistic_end_step) && (steps % RDF_Interval_Steps == 0))
	{
		int i, j;
		double k;
		for (i = 0;i < Particle_Number; i++)
			for (j = i + 1; j < Particle_Number; j++)
			{
				sta_rdf [(int)(Distance_3D_1 (p_head[i].Position[0], p_head[i].Position[1], p_head[i].Position[2], p_head[j].Position[0], p_head[j].Position[1], p_head[j].Position[2], bx, by, bz) * 100.0)] += 1.0;
				//printf("%d\n", (int)(Distance_3D_1 (p_head[i].Position[0], p_head[i].Position[1], p_head[i].Position[2], p_head[j].Position[0], p_head[j].Position[1], p_head[j].Position[2], bx, by, bz) * 100.0));
			}


		if (steps == Statistic_end_step)
		{
			char str_head[80] = "Result_RDF_Output";
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
            FILE *fp_rdf;
            fp_rdf = fopen (str_head, "w");
            for (k = 0.0; k < (float)(RDF_Sta_Length); k += 1.0)
            {
            	//fprintf (fp_rdf, "%lf %lf\n", k * 0.01, sta_rdf[(int)(k)]/(Statistic_end_step - Statistic_start_step + 1.0)/(RDF_Interval_Steps)/(4.0 * Pi * k * 0.01 * k * 0.01 * 0.01)/(0.5 * Particle_Number * (Particle_Number*1.0 - 1.0)));
            	fprintf (fp_rdf, "%lf %lf\n", k * 0.01, sta_rdf[(int)(k)]/(4.0 * Pi * k * 0.01 * k * 0.01 * 0.01)/(0.5 * Particle_Number * (Particle_Number*1.0 - 1.0))*(RDF_Interval_Steps*1.0)/(Statistic_end_step*1.0 - Statistic_start_step*1.0 + 1.0));
            }
            fclose (fp_rdf);
		}
	}
}

void Output_Specific_Configuration_To_File_Mpi (int id, int specific_step, int step, Particle_Model *p_head, double bx, double by, double bz)
{
    if (specific_step == step)
    {
        int total_number, i;
        total_number = Particle_Number;
        char str_head[150] = "Configuration";
        char str_1[20];
        char mpiid[15];
        sprintf (mpiid, "%d", id);
        sprintf (str_1, "%d", specific_step);
        strcat (str_head, "_");
        strcat (str_head, str_1);
        strcat (str_head, "_Mpi_id_");
        strcat (str_head, mpiid); 
        strcat (str_head, ".gro");
        FILE *fp;
        fp = fopen (str_head, "w");
        fprintf (fp, "This is a MC file!\n");
        fprintf (fp, "%d\n", total_number);
        for (i = 0; i < total_number; i++)
            fprintf(fp,"%5d%5s%5s%5d%8.3f%8.3f%8.3f\n", i+1, "SOL", "A", i+1, (p_head + i) -> Position[0], (p_head + i) -> Position[1], (p_head + i) -> Position[2]);

        fprintf (fp, "%8.3f %8.3f %8.3f\n", bx, by, bz);
        fclose (fp);
    }    
}

void Output_Terminal_Configuration_To_File (Particle_Model *p_head, double bx, double by, double bz)
{
    int total_number, i;
    total_number = Particle_Number;
    char str_head[80] = "Result_Terminal_Configuration";
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
    strcat (str_head, ".gro"); 

    FILE *fp;
    fp = fopen (str_head, "w");
    fprintf (fp, "This is a MC file!\n");
    fprintf (fp, "%d\n", total_number);
    for (i = 0; i < total_number; i++)
        fprintf(fp,"%5d%5s%5s%5d%8.3f%8.3f%8.3f\n", i+1, "SOL", "A", i+1, (p_head + i) -> Position[0], (p_head + i) -> Position[1], (p_head + i) -> Position[2]);

    fprintf (fp, "%8.3f %8.3f %8.3f\n", bx, by, bz);
    fclose (fp);
}

void Output_Box_Length_To_File_or_Screen (int type, int interval_steps, int steps, double bl)
{
    /** Type 0: to File!   **
     ** Type 1: to Screen! **
     ** Type 2: Without output **
     ** Interval_step: Output the Box_Length per interval steps **
     **/
    if (type == 0)
    {
        if (steps % interval_steps == 0)
        {
            char str_head[80] = "Result_Box_Length_Output";
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
            FILE *fp_bl;
            fp_bl = fopen (str_head, "a");
            fprintf (fp_bl, "%d %lf\n", steps, bl);
            fclose (fp_bl);
        }
        
    }

    if (type == 1)
    {
        if (steps % interval_steps == 0)
            printf ("%d %lf\n",steps, bl);
    }

    if (type == 2)
    {

    }

}

void Output_Energy_To_File_or_Screen_Mpi (int id, int type, int interval_steps, int steps, double energy)
{
    /** Type 0: to File!   **
     ** Type 1: to Screen! **
     ** Type 2: Without Output **
     ** Interval_step: Output the Energy per interval steps **
     **/
    if (type == 0)
    {
        if (steps % interval_steps == 0)
        {
            char str_head[150] = "Result_Total_Enthalpy_Output";
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
            FILE *fp_energy;
            fp_energy = fopen (str_head, "a");
            fprintf (fp_energy, "%d %lf\n", steps, energy);
            fclose (fp_energy);
        }
        
    }

    if (type == 1)
    {
        if (steps % interval_steps == 0)
            printf ("%d %lf %d\n",steps, energy, id);
    }
    if (type == 2)
    {
        //Do nothing
        
    }
}

void Output_Box_Length_To_File_or_Screen_Mpi (int id, int type, int interval_steps, int steps, double bl)
{
    /** Type 0: to File!   **
     ** Type 1: to Screen! **
     ** Type 2: Without output **
     ** Interval_step: Output the Box_Length per interval steps **
     **/
    if (type == 0)
    {
        if (steps % interval_steps == 0)
        {
            char str_head[80] = "Result_Box_Length_Output";
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
            FILE *fp_bl;
            fp_bl = fopen (str_head, "a");
            fprintf (fp_bl, "%d %lf\n", steps, bl);
            fclose (fp_bl);
        }
        
    }

    if (type == 1)
    {
        if (steps % interval_steps == 0)
            printf ("%d %lf %d\n",steps, bl, id);
    }

    if (type == 2)
    {
        //Do nothing for the test!!
    }

}

#if 0
void Output_Xtcfile_To_File_Or_Screen_Mpi (int id, int interval_steps, int steps, Particle_Model *p_head, XDRFILE *xtc, int natoms, float prec, float tm, double bx, double by, double bz)
{
    if (steps % interval_steps == 0)
    {
        matrix box;
        rvec *x;
        box[0][0] = (float)(bx); box[1][1] = (float)(by); box[2][2] = (float)(bz);
        x = calloc(natoms, sizeof(x[0]));
        int i, j;
        for (i = 0; i < natoms; i++)
            for (j = 0; j < 3; j++)
                x[i][j] = (float)((p_head + i) -> Position[j]);
        write_xtc(xtc, natoms, steps, tm, box, &(x[0]), prec);
        free(x);
    }
}
#endif



#endif
