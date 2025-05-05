#!/bin/bash                                                                                              
if [ "$#" -ne 6 ]; then                                                                                  
    echo "Illegal number of parameters"                                                                  
    exit                                                                                                 
fi                                                                                                       
###Remove dots for name of run
SIGTHname=$(echo "$2"| sed -e 's/\.//g')
Bname=$(echo "$3"| sed -e 's/\.//g')
LAMBTHname=$(echo "$4"| sed -e 's/\.//g')

###Shorten name of damage criterion in name of run
if [ "$6" = "hayurst" ]; then
   CRITname="Ha"
elif [ "$6" = "von mises" ]; then
   CRITname="vM"
elif [ "$6" = "coulomb" ]; then
   CRITname="Co"
elif [ "$6" = "max principal stress" ]; then
   CRITname="MPS"
fi

# ECHO the job that is going to be launch
echo "Launch job: MyTeterousse_Prono_Damage_CavityOnly_$5_B${Bname}_Sigth${SIGTHname}_Lambdah${LAMBTHname}_Crit${CRITname}"                                                                             
# SIF FILE                                                                                               
sed  "s/<SIGTHname>/${SIGTHname}/g;s/<Bname>/${Bname}/g;s/<LAMBTHname>/${LAMBTHname}/g;s/<CRITname>/${CRITname}/g;s/<SIGTH>/$2/g;s/<B>/$3/g;s/<LAMBTH>/$4/g;s/<CRIT>/$6/g" ../SIFs/MyTeterousse_Prono_Damage_CavityOnly_$5.sif  > MyTeterousse_Prono_Damage_CavityOnly_$5_Combi$1_Crit${CRITname}.sif                                                                             
# OAR FILE                                                                                             
sed  "s/<NAME>/MyTeterousse_Prono_Damage_CavityOnly_$5_Combi$1_Crit${CRITname}/g" ../script.oar > script_CavityOnly_$5_Combi$1_Crit${CRITname}.oar                                                                                           

chmod u+x script_CavityOnly_$5_Combi$1_Crit${CRITname}.oar
#### run job                                                                                                 
oarsub -S ./script_CavityOnly_$5_Combi$1_Crit${CRITname}.oar   
