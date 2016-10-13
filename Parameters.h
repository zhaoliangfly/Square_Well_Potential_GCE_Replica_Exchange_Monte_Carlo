#define Coordinate_Number 3
#define BETA_RT 1.0 

/* R = 8.3114
 * 300K RT = 2.493 (RT-1: 0.4011) 
 * 200K RT = 1.662 (RT-1: 0.6017)
 * 100K RT = 0.831 (RT-1: 1.2033)
 */

#define Pi 3.1415926

/* Parameters for Interactive Potential **/
#define Epsilon  1.0  /* KJ/mol */
#define Sigma_Sphere  0.5 /* radius, unit--nm */
#define F_Sigma_Sphere 1.5

#define Initial_Position_Max 10000
#define Particle_Number 300
#define Box_X 7.7
#define Box_Y 7.7
#define Box_Z 7.7

/* Parameters for the Output Control */
#define Statistic_start_step 0
#define Statistic_end_step 5000

/** During prequilibrium, the acceptance rate will be adjusted **/
#define Prequilibrium_Total_Steps 0

#define P_Reference 6e-02

#define P_bias 0.95
#define N_Volume_Steps 80

//#define N_AVB_Steps 100


#define V_max  0.05


/** GCE_Common_Single_Particle_Metropolis_Sampling **/
#define Max_Seeding_Time 1000  /** int type **/
#define Total_Steps 300000000 /** int type **/
#define Output_Interval_Steps 1000
#define RDF_Interval_Steps 10000
#define RDF_Sta_Length 15000
#define Adjust_Ratio_Interval_Steps 1000
#define Object_Ratio 0.5
#define Delta_R 0.05  //distance for the translation, unit nm. 

/** Parameters for GCE **/
//#define H0 0.0 
#define Alpha 6e-03 

/** Output to file or screen control **
 ** Type 0: to File!                 **
 ** Type 1: to Screen!               **
 ** Type 2: Without Output!          **
 ** Interval_step: Output the Energy per interval steps **
 **/
#define To_File_Or_Screen  2
#define To_File_Or_Screen_Mpi_Enthalpy 0
#define To_File_Or_Screen_Mpi_Box_Length 0
#define To_File_Or_Screen_Mpi_Acceptance_Ratio 0
#define To_File_Or_Screen_Mpi_Swap_Acceptance_Ratio 0

/** Replicas interval **/
#define Replicas_Exchange_Interval_Steps 100
#define Calculate_Swap_Acceptance_Ratio_Interval_Steps 1000000

/** Write_Xtcfile **/
#define Write_Xtcfile 0
#define Xtc_Interval_Steps 1000
