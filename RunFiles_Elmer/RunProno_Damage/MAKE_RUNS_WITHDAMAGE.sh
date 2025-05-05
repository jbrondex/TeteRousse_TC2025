#!/bin/bash
cmin=1
cmax=8
step="Refill20122013" ## "Refill20102011" ##"Pumping2010"

##declare list of investigated damage criterion
declare -a Criterion=("hayurst" "von mises" "coulomb" "max principal stress")

## get the length of the array
Criterionlength=${#Criterion[@]}
##Start loop on combinations of damage parameters
for c in $(seq $cmin $cmax); do 
   DIRECTORY=RunProno_Damage_Combi$c
   if [ -d "$DIRECTORY" ]; then
      echo "DIRECTORY ALREADY EXISTING: " $DIRECTORY
      echo "ABORT"
      exit
   fi

#### Create run dir
   mkdir $DIRECTORY
   cd $DIRECTORY
#### Create output dir
   mkdir output
#### Create dir for .dat.names of SaveLine (bug of Solver ?)
   mkdir ScalarOutput
#### copy Data directory
   cp -r ../Data ./
#### copy mesh
   mkdir teterousse
   cp -r ../teterousse/partitioning.8/ teterousse/.
   cp ../teterousse/RunSteadyInit_InitGeoAndGMAndT_.result.* teterousse/.
#### For each damage criterion, create sif and launch run for 3 pressure case
   for (( i=0; i<${Criterionlength}; i++ )); 
   do 
   	sh ../SCRIPTs/Run_Damage_CavityOnly.sh $c $(awk -v n=$c 'NR == n' ../Dam_Sigmath_B_Lambdah.IN) $step "${Criterion[$i]}"
   done
###
   cd ../
done

