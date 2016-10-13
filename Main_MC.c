#include "stdio.h"
#include "math.h"
#include "time.h"
#include "Parameters.h"
#include "mpi.h"
//#include "xdrfile.h"

#define AVB_Algorithm 0

/* Definition of the Structures to Describe the Particles */
// Attention! For the convienience of MPI programming, we change the struct to 2d array!
struct data
{
	double Position[Coordinate_Number];
};

typedef struct data Particle_Model;
Particle_Model Particles[Particle_Number];

//double Particles[Particle_Number][Coordinate_Number];

/** To save the computation resource, the variant is introduced whether the configuration is necessary to be recomputed! **/
double H0;
double H0_Interval = 25.0;
double H0_Base = -1550.0;
double former_energy;

//#if Write_Xtcfile
//int natoms = Particle_Number;
//int magic_number = 1995;
//float prec, tm;
//XDRFILE *xtc;
//#endif

#include "Initialization_System.h"
#include "Square_Well_Potential_Definition.h"
#include "GCE_Metropolis_Move_For_NPT.h"
#include "GCE_Volume_Move_For_NPT.h"
//#include "GCE_AVB_Move_For_NPT.h"
#include "Adjust_Acceptance_Ratio.h"
#include "Output_Energy_Configuration.h"
#include "Swap_Acceptance_Ratio.h"

int main(int argc, char *argv[])
{
    int myid, numprocs;
    int namelen;
    char processor_name[MPI_MAX_PROCESSOR_NAME];
    MPI_Status status;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Get_processor_name(processor_name, &namelen);

    H0 = H0_Interval * (double)(myid) + H0_Base;// Different H0 for a series of simulations. 

	int steps = 0;
	double b_x, b_y, b_z;
	b_x = Box_X; b_y = Box_Y; b_z = Box_Z; //NPT will change the length of box

	Initialize_Particle_Model(1, &(Particles[0]), b_x, b_y, b_z);

	/** Write the first frame to Xtcfile **/
//#if Write_Xtcfile
//	char str_head[150] = "traj";
//    char mpiid[15];
//    sprintf(mpiid, "%d", myid);
//    strcat(str_head, "_");
//    strcat(str_head, mpiid);
//    strcat(str_head, ".xtc");
//   xtc = xdrfile_open(str_head, "w");

//	Output_Xtcfile_To_File_Or_Screen_Mpi(myid, Xtc_Interval_Steps, steps, &(Particles[0]), xtc, natoms, prec, tm, b_x, b_y, b_z);
//#endif

    former_energy = Total_Energy_Square_Well_Potential(&(Particles[0]), b_x, b_y, b_z);

	double Accept_number, Reject_number, Accept_number_swap, Reject_number_swap;
	double delta_r = Delta_R;
	Accept_number = 0.0; Reject_number = 0.0; Accept_number_swap = 0.0; Reject_number_swap = 0.0;

	#if 0
	double rdf[RDF_Sta_Length];
	int l = 0;
	for (l = 0; l < RDF_Sta_Length; l++)
		rdf[l] = 0.0;
	#endif

	dsfmt_t dsfmt;
	//Attention!!!!! In the Mpi computation, the seed of random generator must be different.
	//So, we change the (unsigned)time(NULL) to myid.
	//dsfmt_init_gen_rand (&dsfmt, (unsigned)time(NULL));
	dsfmt_init_gen_rand(&dsfmt, myid + (unsigned)time(NULL));
	
	//while (steps <= Prequilibrium_Total_Steps)

	while (steps <= Total_Steps)
	{
		int rand_number;

#if AVB_Algorithm //AVB has not been finished!
		rand_number = (int)(dsfmt_genrand_open_open(&dsfmt) * (Particle_Number + N_Volume_Number + N_AVB_Number));
		if (rand_number < Particle_Number)  
			GCE_Metropolis_Move_For_NPT_Function(&steps, &(Particles[0]), &dsfmt, &delta_r, &Accept_number, &Reject_number, b_x, b_y, b_z);
		else if (rand_number < Particle_Number + N_Volume_Number)
			GCE_Volume_Move_For_NPT_Function(&steps, &(Particles[0]), &dsfmt, &delta_r, &Accept_number, &Reject_number, &b_x, &b_y, &b_z);
		else
			GCE_AVB_Move_For_NPT_Function(&steps, &(Particles[0]), &dsfmt, &delta_r, &Accept_number, &Reject_number, &b_x, &b_y, &b_z);
#else
		rand_number = (int)(dsfmt_genrand_open_open(&dsfmt) * 90.0) + 10;
		if (rand_number <= 90)   // 90% 
			GCE_Metropolis_Move_For_NPT_Function_Mpi(myid, &steps, &(Particles[0]), &dsfmt, &delta_r, &Accept_number, &Reject_number, b_x, b_y, b_z);
		else // 10%
			GCE_Volume_Move_For_NPT_Function_Mpi(myid, &steps, &(Particles[0]), &dsfmt, &delta_r, &Accept_number, &Reject_number, &b_x, &b_y, &b_z);
#endif
		
		MPI_Barrier(MPI_COMM_WORLD);

		/** Codes below are specially for the replica-exchange. **
		 **                                                     **
		 **                                                     **
		 **/

		//if (steps == Replicas_Exchange_Steps)
		if (steps % Replicas_Exchange_Interval_Steps == 0)
		{
			double Random_For_Replica_Methods;
			if (myid == 0)
				Random_For_Replica_Methods = dsfmt_genrand_open_open(&dsfmt);
			MPI_Bcast(&Random_For_Replica_Methods, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
			
			if (Random_For_Replica_Methods > 0.5) //Method 1
			{
				//===========Method 1: (0, 1) (2, 3) (4, 5)......==========//
				double random_exchange;
				random_exchange = dsfmt_genrand_open_open(&dsfmt);
				if (myid % 2 == 0) //We must ensure the total number of processes is even.
					MPI_Send(&(random_exchange), 1, MPI_DOUBLE, myid + 1, 99, MPI_COMM_WORLD);
				else
					MPI_Recv(&(random_exchange), 1, MPI_DOUBLE, myid - 1, 99, MPI_COMM_WORLD, &status);
				MPI_Barrier(MPI_COMM_WORLD);

				double H, E, temp_H;
				H = former_energy + P_Reference * b_x * b_y * b_z;
				E = former_energy; //Attention!! Don't forget to exchange the E, as the E is calculated on the basis of former energy.

				if (myid % 2 == 1)
					MPI_Send(&H, 1, MPI_DOUBLE, myid - 1, 98, MPI_COMM_WORLD);
				else
					MPI_Recv(&temp_H, 1, MPI_DOUBLE, myid + 1, 98, MPI_COMM_WORLD, &status);
				MPI_Barrier(MPI_COMM_WORLD);
			
				if (myid % 2 == 0)
					H = temp_H - H;

				if (myid % 2 == 0)//delta_H is stored in H.
					MPI_Send(&H, 1, MPI_DOUBLE, myid + 1, 97, MPI_COMM_WORLD);
				else
					MPI_Recv(&H, 1, MPI_DOUBLE, myid - 1, 97, MPI_COMM_WORLD, &status);
				MPI_Barrier(MPI_COMM_WORLD);

				double temp_position[Particle_Number][Coordinate_Number];
				int i, j;
				for (i = 0; i < Particle_Number; i++)
				{
					for (j = 0; j < Coordinate_Number; j++)
						temp_position[i][j] = Particles[i].Position[j];
				}
				double temp_box[3];
				temp_box[0] = b_x; temp_box[1] = b_y; temp_box[2] = b_z;
				//printf("Steps: %d %d\n", steps, myid);

				//Begin to exchange the E.
				if (myid % 2 == 0 && random_exchange < exp(-1.0 * Alpha * H * H0_Interval))
					MPI_Send(&E, 1, MPI_DOUBLE, myid + 1, 80, MPI_COMM_WORLD); //Send: From Even to Odd
				if (myid % 2 == 1 && random_exchange < exp(-1.0 * Alpha * H * H0_Interval))
					MPI_Recv(&former_energy, 1, MPI_DOUBLE, myid - 1, 80, MPI_COMM_WORLD, &status); //Receive: From Even to Odd
				MPI_Barrier(MPI_COMM_WORLD);

				if (myid % 2 == 1 && random_exchange < exp(-1.0 * Alpha * H * H0_Interval))
					MPI_Send(&E, 1, MPI_DOUBLE, myid - 1, 79, MPI_COMM_WORLD); //Send: From Even to Odd
				if (myid % 2 == 0 && random_exchange < exp(-1.0 * Alpha * H * H0_Interval))
					MPI_Recv(&former_energy, 1, MPI_DOUBLE, myid + 1, 79, MPI_COMM_WORLD, &status); //Receive: From Even to Odd
				MPI_Barrier(MPI_COMM_WORLD);


				//Begin to exchange the coordinates and box.
				//Pay attentions to the Swap Acceptance Ratio!
				if (myid % 2 == 0 && random_exchange < exp(-1.0 * Alpha * H * H0_Interval))
					Accept_number_swap += 1.0;
				if (myid % 2 == 0 && random_exchange > exp(-1.0 * Alpha * H * H0_Interval))
					Reject_number_swap += 1.0;

				if (myid % 2 == 0 && random_exchange < exp(-1.0 * Alpha * H * H0_Interval))//Send: From Even to Odd
				{
					for (i = 0; i < Particle_Number; i++)
						MPI_Send(&(temp_position[i][0]), 3, MPI_DOUBLE, myid + 1, i, MPI_COMM_WORLD);
					MPI_Send(&(temp_box[0]), 1, MPI_DOUBLE, myid + 1, 96, MPI_COMM_WORLD);
					MPI_Send(&(temp_box[1]), 1, MPI_DOUBLE, myid + 1, 95, MPI_COMM_WORLD);
					MPI_Send(&(temp_box[2]), 1, MPI_DOUBLE, myid + 1, 94, MPI_COMM_WORLD);
				}
				if (myid % 2 == 1 && random_exchange < exp(-1.0 * Alpha * H * H0_Interval)) //Receive: From Even to Odd
				{
					for (i = 0; i < Particle_Number; i++)
						MPI_Recv(&(Particles[i].Position[0]), 3, MPI_DOUBLE, myid - 1, i, MPI_COMM_WORLD, &status);
					MPI_Recv(&b_x, 1, MPI_DOUBLE, myid - 1, 96, MPI_COMM_WORLD, &status);
					MPI_Recv(&b_y, 1, MPI_DOUBLE, myid - 1, 95, MPI_COMM_WORLD, &status);
					MPI_Recv(&b_z, 1, MPI_DOUBLE, myid - 1, 94, MPI_COMM_WORLD, &status);
				}
				MPI_Barrier(MPI_COMM_WORLD);

				if (myid % 2 == 1 && random_exchange < exp(-1.0 * Alpha * H * H0_Interval)) //Send: From Odd to Even 
				{
					for (i = 0; i < Particle_Number; i++)
						MPI_Send(&(temp_position[i][0]), 3, MPI_DOUBLE, myid - 1, i, MPI_COMM_WORLD);
					MPI_Send(&(temp_box[0]), 1, MPI_DOUBLE, myid - 1, 93, MPI_COMM_WORLD);
					MPI_Send(&(temp_box[1]), 1, MPI_DOUBLE, myid - 1, 92, MPI_COMM_WORLD);
					MPI_Send(&(temp_box[2]), 1, MPI_DOUBLE, myid - 1, 91, MPI_COMM_WORLD);
				}
				if (myid % 2 == 0 && random_exchange < exp(-1.0 * Alpha * H * H0_Interval)) //Receive: From Odd to Even
				{
					for (i = 0; i < Particle_Number; i++)
						MPI_Recv(&(Particles[i].Position[0]), 3, MPI_DOUBLE, myid + 1, i, MPI_COMM_WORLD, &status);
					MPI_Recv(&b_x, 1, MPI_DOUBLE, myid + 1, 93, MPI_COMM_WORLD, &status);
					MPI_Recv(&b_y, 1, MPI_DOUBLE, myid + 1, 92, MPI_COMM_WORLD, &status);
					MPI_Recv(&b_z, 1, MPI_DOUBLE, myid + 1, 91, MPI_COMM_WORLD, &status);
				}
				MPI_Barrier(MPI_COMM_WORLD);
			}
			else
			{
				//===========Method 2: 0, (1, 2), (3, 4), 5......==========//
				double random_exchange;
				random_exchange = dsfmt_genrand_open_open(&dsfmt);
				if (myid % 2 == 1 && myid < numprocs - 1) //We must ensure the total number of processes is even.
					MPI_Send(&(random_exchange), 1, MPI_DOUBLE, myid + 1, 99, MPI_COMM_WORLD);
				if (myid % 2 == 0 && myid != 0)
					MPI_Recv(&(random_exchange), 1, MPI_DOUBLE, myid - 1, 99, MPI_COMM_WORLD, &status);
				MPI_Barrier(MPI_COMM_WORLD);

				double H, E, temp_H;
				H = former_energy + P_Reference * b_x * b_y * b_z;
				E = former_energy; //Attention!! Don't forget to exchange the E, as the E is calculated on the basis of former energy.
				
				if (myid % 2 == 0 && myid != 0)
					MPI_Send(&H, 1, MPI_DOUBLE, myid - 1, 98, MPI_COMM_WORLD);
				if (myid % 2 == 1 && myid < numprocs - 1)
					MPI_Recv(&temp_H, 1, MPI_DOUBLE, myid + 1, 98, MPI_COMM_WORLD, &status);
				MPI_Barrier(MPI_COMM_WORLD);

				if (myid % 2 == 1 && myid < numprocs - 1)
					H = temp_H - H;

				if (myid % 2 == 1 && myid < numprocs - 1)//delta_H is stored in H.
					MPI_Send(&H, 1, MPI_DOUBLE, myid + 1, 97, MPI_COMM_WORLD);
				if (myid % 2 == 0 && myid != 0)
					MPI_Recv(&H, 1, MPI_DOUBLE, myid - 1, 97, MPI_COMM_WORLD, &status);
				MPI_Barrier(MPI_COMM_WORLD);

				double temp_position[Particle_Number][Coordinate_Number];
				int i, j;
				for (i = 0; i < Particle_Number; i++)
				{
					for (j = 0; j < Coordinate_Number; j++)
						temp_position[i][j] = Particles[i].Position[j];
				}
				double temp_box[3];
				temp_box[0] = b_x; temp_box[1] = b_y; temp_box[2] = b_z;
				//printf("Steps: %d %d\n", steps, myid);

				//Begin to exchange the E.
				if (myid % 2 == 1 && myid < numprocs - 1 && random_exchange < exp(-1.0 * Alpha * H * H0_Interval))
					MPI_Send(&E, 1, MPI_DOUBLE, myid + 1, 80, MPI_COMM_WORLD); //Send: From Odd to Even
				if (myid % 2 == 0 && myid != 0 && random_exchange < exp(-1.0 * Alpha * H * H0_Interval))
					MPI_Recv(&former_energy, 1, MPI_DOUBLE, myid - 1, 80, MPI_COMM_WORLD, &status); //Receive: From Odd to Even
				MPI_Barrier(MPI_COMM_WORLD);

				if (myid % 2 == 0 && myid != 0 && random_exchange < exp(-1.0 * Alpha * H * H0_Interval))
					MPI_Send(&E, 1, MPI_DOUBLE, myid - 1, 79, MPI_COMM_WORLD); //Send: From Even to Odd
				if (myid % 2 == 1 && myid < numprocs - 1 && random_exchange < exp(-1.0 * Alpha * H * H0_Interval))
					MPI_Recv(&former_energy, 1, MPI_DOUBLE, myid + 1, 79, MPI_COMM_WORLD, &status); //Receive: From Even to Odd
				MPI_Barrier(MPI_COMM_WORLD);

				//Begin to exchange the coordinates and box.
				if (myid % 2 == 1 && myid < numprocs - 1 && random_exchange < exp(-1.0 * Alpha * H * H0_Interval))
					Accept_number_swap += 1.0;
				if (myid % 2 == 1 && myid < numprocs - 1 && random_exchange > exp(-1.0 * Alpha * H * H0_Interval))
					Reject_number_swap += 1.0;

				if (myid % 2 == 1 && myid < numprocs - 1 && random_exchange < exp(-1.0 * Alpha * H * H0_Interval))//Send: From Odd to Even
				{
					for (i = 0; i < Particle_Number; i++)
						MPI_Send(&(temp_position[i][0]), 3, MPI_DOUBLE, myid + 1, i, MPI_COMM_WORLD);
					MPI_Send(&(temp_box[0]), 1, MPI_DOUBLE, myid + 1, 96, MPI_COMM_WORLD);
					MPI_Send(&(temp_box[1]), 1, MPI_DOUBLE, myid + 1, 95, MPI_COMM_WORLD);
					MPI_Send(&(temp_box[2]), 1, MPI_DOUBLE, myid + 1, 94, MPI_COMM_WORLD);
				}
				if (myid % 2 == 0 && myid != 0 && random_exchange < exp(-1.0 * Alpha * H * H0_Interval)) //Receive: From Odd to Even
				{
					for (i = 0; i < Particle_Number; i++)
						MPI_Recv(&(Particles[i].Position[0]), 3, MPI_DOUBLE, myid - 1, i, MPI_COMM_WORLD, &status);
					MPI_Recv(&b_x, 1, MPI_DOUBLE, myid - 1, 96, MPI_COMM_WORLD, &status);
					MPI_Recv(&b_y, 1, MPI_DOUBLE, myid - 1, 95, MPI_COMM_WORLD, &status);
					MPI_Recv(&b_z, 1, MPI_DOUBLE, myid - 1, 94, MPI_COMM_WORLD, &status);
				}
				MPI_Barrier(MPI_COMM_WORLD);

				if (myid % 2 == 0 && myid != 0 && random_exchange < exp(-1.0 * Alpha * H * H0_Interval)) //Send: From Even to Odd 
				{
					for (i = 0; i < Particle_Number; i++)
						MPI_Send(&(temp_position[i][0]), 3, MPI_DOUBLE, myid - 1, i, MPI_COMM_WORLD);
					MPI_Send(&(temp_box[0]), 1, MPI_DOUBLE, myid - 1, 93, MPI_COMM_WORLD);
					MPI_Send(&(temp_box[1]), 1, MPI_DOUBLE, myid - 1, 92, MPI_COMM_WORLD);
					MPI_Send(&(temp_box[2]), 1, MPI_DOUBLE, myid - 1, 91, MPI_COMM_WORLD);
				}
				if (myid % 2 == 1 && myid < numprocs - 1 && random_exchange < exp(-1.0 * Alpha * H * H0_Interval)) //Receive: From Even to Odd
				{
					for (i = 0; i < Particle_Number; i++)
						MPI_Recv(&(Particles[i].Position[0]), 3, MPI_DOUBLE, myid + 1, i, MPI_COMM_WORLD, &status);
					MPI_Recv(&b_x, 1, MPI_DOUBLE, myid + 1, 93, MPI_COMM_WORLD, &status);
					MPI_Recv(&b_y, 1, MPI_DOUBLE, myid + 1, 92, MPI_COMM_WORLD, &status);
					MPI_Recv(&b_z, 1, MPI_DOUBLE, myid + 1, 91, MPI_COMM_WORLD, &status);
				}
				MPI_Barrier(MPI_COMM_WORLD);
			}

		}//End the steps == Replicas_Exchange_Steps

        MPI_Barrier(MPI_COMM_WORLD);
        Output_Swap_Acceptance_Ratio_Mpi(myid, numprocs, To_File_Or_Screen_Mpi_Swap_Acceptance_Ratio, steps, &Accept_number_swap, &Reject_number_swap);	

//#if Write_Xtcfile
//		Output_Xtcfile_To_File_Or_Screen_Mpi(myid, Xtc_Interval_Steps, steps, &(Particles[0]), xtc, natoms, prec, tm, b_x, b_y, b_z);	
//#endif
	
		//Output_Radius_Distribution_Function_To_File (steps, &(rdf[0]), &(Particles[0]), b_x, b_y, b_z);
		Adjust_Acceptance_Mpi(myid, To_File_Or_Screen_Mpi_Acceptance_Ratio, Adjust_Ratio_Interval_Steps, steps, &Accept_number, &Reject_number, &delta_r);
	    Output_Box_Length_To_File_or_Screen_Mpi(myid, To_File_Or_Screen_Mpi_Box_Length, Output_Interval_Steps, steps, b_x);

	    //Output_Specific_Configuration_To_File_Mpi(myid, 1000, steps, &(Particles[0]), b_x, b_y, b_z);
	}
    
    //Output_Terminal_Configuration_To_File (&(Particles[0]), b_x, b_y, b_z);

//#if Write_Xtcfile
//    xdrfile_close(xtc);
//#endif

    MPI_Finalize();

	return 0;
}
