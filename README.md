# How to model crevasse initiation ? Lessons from the artificial drainage of a water-filled cavity on the Tête Rousse Glacier (Mont Blanc range, France)

**Authors:** Julien Brondex; Olivier Gagliardini, Adrien Gilbert, and Emmanuel Thibert

**Contact author:** Julien Brondex - julien.brondex@univ-grenoble-alpes.fr 

*Univ. Grenoble Alpes, CNRS, IRD, Grenoble INP, INRAE, IGE, 38000 Grenoble, France*

## Motivations

Our primary goal is to establish a correlation between the occurrence of crevasses in glaciers and the distribution and magnitude of stresses. To this end, we use the carefully monitored artificial drainage of a water-filled cavity on Tête Rousse Glacier (Mont Blanc range, France) to constrain the failure criterion and stress threshold that best reproduce the field of circular crevasses mapped at the surface during the summer following the pumping operation. We performed numerical simulations using the finite element code Elmer/Ice. Since the evolution of the water level and the cavity volume (and, to a lesser extent, its geometry) are well-constrained, water pressure against the cavity wall can be prescribed as a time-dependent boundary condition, and stress fields are inferred directly from the force balance. In other words, our numerical experiments are entirely force-driven. The resulting deformations serve as independent control variables to assess the nature of the mechanical response. This work is at the origin of an article submitted to The Cryosphere on 2025/05/07 and accepted on 2026/03/06 

## Quick install
For numerical simulations, we used version 9.0 of Elmer/Ice (Rev:0f18c8f). The Elmer/Ice code is publicly available through GitHub (https://github.com/ElmerCSC/elmerfem last access: Mai 2025). 

For postprocessing, we used the following Python packages:

python==3.12

numpy==2.2.4

matplotlib==3.10.1

pandas==2.2.3

scipy==1.15.2

## Structure github repository

```
├── README.md                           <- File describing the github repository. You are currently reading this file.
│
├── Data                                <- Contains all field data used in the study                            
│   ├── Crevasses                       <- Contains file with positions of crevasses
│   :    
│   ├── SMB                             <- Contains file of reconstructed daily SMB             
│   :                                                
│   ├── Stakes                          <- Contains files with positions of stakes used to monitor the displacement of the cavity roof in 2010, 2011, 2012
│   :                                                
│   ├── Temperature                     <- Contains raw data of thermistor measurements and python script to produce the reconstructed temperature field 
│   :                                                
│   ├── Topo                            <- Contains all files required to construct the initial geometry of the computational domains
│   :                                                
│   ├── WaterLevel                      <- Contains raw data as well as reconstructed water level within the cavity
│
│
├── PostProcessing                      <- Contains all python scripts used to produce the figures of the paper
│
├── RunFiles_Elmer                      <- Contains all files required to run the numerical simulations, provided that a compatible version of Elmer/Ice is installed               
│   ├── bin                             <- An empty folder in which binaries resulting from the compilation of fortran files in SRC through Makefile will be stored 
│   :                                
│   ├── RunDiag                         <- All input files to run the diagnostic simulations presented in the paper (Elastic, linear viscous, non linear Glen) 
│   :                                
│   ├── Run_Initialisation              <- All files needed to generate the initial mesh and run a simulation that initializes the geometry and some fields (e.g., temperature)  
│   :                                
│   ├── RunProno_Damage                 <- All files to run all simulations with the 8 combinations of damage parameters with the 4 tested criteria automatically  
│   :                                
│   ├── RunProno_NoDamage               <- All files to run pronostic simu without damage with sensitivity tests on temperature field and the way water pressure is applied  
│   :                                
│   ├── SRC                             <- Source code that is not available in the distribution but was developped specifically for this study 
│   :                                
```
## Run the simulations

The following provides details on how to run the simulations presented in this study. It is not intended to serve as a user manual for Elmer. Documentation regarding Elmer can be found on the Elmer git : https://github.com/ElmerCSC/elmerfem (last accessed: March 2026).

All files required to run the simulations are in the RunFiles_Elmer. Of course, all path contained in these files need to be adapted. 

First, compile the source code developped specifically for this study using the Makefile executable. All binaries should be automatically sent in the bin folder.
Second, go in the Run_Initialisation folder and build the mesh using the MakeMSH file. Then, initialise the geometry and fields by running the MyTeterousse_Init_WithT.sif file. 
Third, you can navigate to any of the RunDiag, RunProno_Damage, RunProno_NoDamage folders and run the simulations you are interested in. Note that the mesh folder and result files generated during the initialisation run must be copied into the folder where you intend to run the simulation. Note also that the RunProno_Damage folder is structured to automatically run simulations for the 8 combinations of damage parameters tested and the 4 damage criteria evaluated. This is achieved by executing the MAKE_RUNS_WITHDAMAGE.sh script. The simulation step you wish to run ("Pumping2010", "Refill20102011", "Pumping2011", and so on) must be specified within this script. The script is configured for use with the OAR task manager and should be adapted to your own job scheduling environment as needed.

