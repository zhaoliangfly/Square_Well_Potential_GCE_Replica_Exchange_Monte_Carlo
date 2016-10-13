#!/bin/bash
PR="0.060"
KT="1.000"
ALPHA="0.006"

#Parameters for Sampling
H0_Base="-1550.000"
H0_Interval="25.0"
H_Translation="5000.0" #Ensure the value of H is positive.
H_Upper="10000.0"
Box_Upper="100000.0"


#Parameters for Statistic
Time_Ini="200000000.0"
Time_Fin="300000000.0"

#The file name
Folders=([0]=Acceptance_Ratio [1]=Box_Length [2]=Swap_Acceptance_Ratio [3]=Total_Enthalpy)
if [ -d "Sta_Enthalpy" ];then
rm -rf Sta_Enthalpy
fi

if [ -d "Sta_Swap_Acceptance_Ratio" ];then
rm -rf Sta_Swap_Acceptance_Ratio
fi

if [ -f "GCE_Temperature_vs_Density_AL_${ALPHA}_BE_${KT}_Pr_${PR}_Interval_25.dat" ];then
rm GCE_Temperature_vs_Density_AL_${ALPHA}_BE_${KT}_Pr_${PR}.dat
fi

if [ -f "GCE_Beta_vs_Enthalpy_AL_${ALPHA}_BE_${KT}_Pr_${PR}_Interval_25.dat" ];then
rm GCE_Beta_vs_Enthalpy_AL_${ALPHA}_BE_${KT}_Pr_${PR}_Interval_25.dat
fi

echo
echo ===== Beta_${KT}, Pressure_${PR}, Alpha_${ALPHA} =====
echo
echo
echo "Starting Box_Length and Enthalpy Calculations..."
echo

mkdir Sta_Enthalpy

#The first calculation: beta(H) with errorbars
for id in {0..71};do

cd ${Folders[3]}
H0_file=`echo "scale=3;${H0_Base}+${H0_Interval}*${id}" | bc`
if [ ${H0_file:0:1} == "-" ]; then
	H0_file=${H0_file:1}
fi
if [ "$H0_file" == "0" ]; then
	H0_file="0.000"
fi
File_name_1="Result_${Folders[3]}_Output_AL_${ALPHA}_BE_${KT}_Pr_${PR}_H0_${H0_file}_Mpi_id_${id}.dat"

H0=`echo "scale=3;${H0_Base}+${H0_Interval}*${id}" | bc`
#average of Enthalpy
Average_H=`awk -v ini="$Time_Ini" -v fin="$Time_Fin"  'BEGIN{Average=0.0;Count=0.0;}{if($1>ini && $1<fin) {Average+=$2;Count+=1.0;}} END{printf("%lf",Average/Count)}' ${File_name_1}`
#variation of Enthalpy
Var_H=`awk -v av_h="$Average_H" -v ini="$Time_Ini" -v fin="$Time_Fin" 'BEGIN{Var=0.0;Count=0.0;}{if($1>ini && $1<fin){Var+=($2-av_h)*($2-av_h);Count+=1.0;}}END{printf("%lf",sqrt(Var/(Count-1.0)))}' ${File_name_1}`

#statistic of Enthalpy
awk -v ini="$Time_Ini" -v fin="$Time_Fin" -v trans="$H_Translation" -v max="$H_Upper" 'BEGIN{count=0.0;for(i=0;i<max;i++) Sta_H[i]=0.0;}{if($1>ini && $1<fin){count+=1.0;Sta_H[int($2+trans)]+=1.0;}} END{for(i=0;i<max;i++){printf("%lf %lf\n", i*1.0-trans, Sta_H[i]/(count-1.0));}}' ${File_name_1} > ./../Sta_Enthalpy/Sta_${Folders[3]}_Output_AL_${ALPHA}_BE_${KT}_Pr_${PR}_H0_${H0_file}_Mpi_id_${id}.dat 

#Beta_value
awk -v h0="${H0}" -v alpha="${ALPHA}" -v beta_0="${KT}" '{printf("%d %lf\n", $1, beta_0+alpha*($2-h0))}' ${File_name_1} > temp_beta.dat
#average of Beta
Average_Beta=`awk -v ini="$Time_Ini" -v fin="$Time_Fin"  'BEGIN{Average=0.0;Count=0.0;}{if($1>ini && $1<fin) {Average+=$2;Count+=1.0;}} END{printf("%lf",Average/Count)}' temp_beta.dat`
#variation of Beta
Var_Beta=`awk -v av_b="$Average_Beta" -v ini="$Time_Ini" -v fin="$Time_Fin" 'BEGIN{Var=0.0;Count=0.0;}{if($1>ini && $1<fin){Var+=($2-av_b)*($2-av_b);Count+=1.0;}}END{printf("%lf",sqrt(Var/(Count-1.0)))}' temp_beta.dat`

rm temp_beta.dat
echo $Average_H $Average_Beta $Var_H $Var_Beta >> ./../GCE_Beta_vs_Enthalpy_AL_${ALPHA}_BE_${KT}_Pr_${PR}_Interval_25.dat

#T_value
awk -v h0="${H0}" -v alpha="${ALPHA}" -v beta_0="${KT}" '{printf("%d %lf\n", $1, 1.0/(beta_0+alpha*($2-h0)))}' ${File_name_1} > temp_t.dat
#average of T
Average_T=`awk -v ini="$Time_Ini" -v fin="$Time_Fin"  'BEGIN{Average=0.0;Count=0.0;}{if($1>ini && $1<fin) {Average+=$2;Count+=1.0;}} END{printf("%lf",Average/Count)}' temp_t.dat`
#variation of T
Var_T=`awk -v av_t="$Average_T" -v ini="$Time_Ini" -v fin="$Time_Fin" 'BEGIN{Var=0.0;Count=0.0;}{if($1>ini && $1<fin){Var+=($2-av_t)*($2-av_t);Count+=1.0;}}END{printf("%lf",sqrt(Var/(Count-1.0)))}' temp_t.dat`

rm temp_t.dat

cd ..

cd ${Folders[1]}
File_name_2="Result_${Folders[1]}_Output_AL_${ALPHA}_BE_${KT}_Pr_${PR}_H0_${H0_file}_Mpi_id_${id}.dat"
#Density_value
awk '{printf("%d %lf\n", $1, 300.0/($2*$2*$2))}' ${File_name_2} > temp_density.dat
#average of Density
Average_Density=`awk -v ini="$Time_Ini" -v fin="$Time_Fin"  'BEGIN{Average=0.0;Count=0.0;}{if($1>ini && $1<fin) {Average+=$2;Count+=1.0;}} END{printf("%lf",Average/Count)}' temp_density.dat`
#variation of Density
Var_Density=`awk -v av_d="$Average_Density" -v ini="$Time_Ini" -v fin="$Time_Fin" 'BEGIN{Var=0.0;Count=0.0;}{if($1>ini && $1<fin){Var+=($2-av_d)*($2-av_d);Count+=1.0;}}END{printf("%lf",sqrt(Var/(Count-1.0)))}' temp_density.dat`

rm temp_density.dat
echo $Average_Density $Average_T $Var_Density $Var_T >> ./../GCE_Temperature_vs_Density_AL_${ALPHA}_BE_${KT}_Pr_${PR}_Interval_25.dat
cd ..

done

mkdir Sta_Swap_Acceptance_Ratio
for id in {0..70};do
cd ${Folders[2]}
H0_file=`echo "scale=3;${H0_Base}+${H0_Interval}*${id}" | bc`
if [ ${H0_file:0:1} == "-" ]; then
	H0_file=${H0_file:1}
fi
if [ "$H0_file" == "0" ]; then
	H0_file="0.000"
fi
File_name_3="Result_${Folders[2]}_AL_${ALPHA}_BE_${KT}_Pr_${PR}_H0_${H0_file}_Mpi_id_${id}.dat"
#average of Swap
Average_Swap=`awk -v ini="$Time_Ini" -v fin="$Time_Fin"  'BEGIN{Average=0.0;Count=0.0;}{if($1>ini && $1<fin) {Average+=$2;Count+=1.0;}} END{printf("%lf",Average/Count)}' ${File_name_3}`
#variation of Swap
Var_Swap=`awk -v av_s="$Average_Swap" -v ini="$Time_Ini" -v fin="$Time_Fin" 'BEGIN{Var=0.0;Count=0.0;}{if($1>ini && $1<fin){Var+=($2-av_s)*($2-av_s);Count+=1.0;}}END{printf("%lf",sqrt(Var/(Count-1.0)))}' ${File_name_3}`

echo $id $Average_Swap $Var_Swap >> ./../Sta_Swap_Acceptance_Ratio/Sta_${Folders[2]}_AL_${ALPHA}_BE_${KT}_Pr_${PR}_Interval_25.dat
cd ..
done


echo
echo "End Box_Length and Enthalpy Calculations..."
echo ===== End =====

