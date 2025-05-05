"""
@author: jbrondex

Description:
------------
This file aims at plotting the evolution of the mean vertical velocity above the cavity between first pumping of 2010 to 2014. We test the following
damage criterion : Max Principle Stress, von Mises, Hayurst, Coulomb

"""

################################################################################
# import #####
################################################################################
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cmx

from pathlib import Path
from datetime import datetime, date, timedelta ###To handle dates
from matplotlib.dates import MonthLocator, DayLocator, HourLocator, DateFormatter, drange
import numpy as np
from matplotlib.ticker import FormatStrFormatter
from matplotlib.dates import date2num
import matplotlib.gridspec as gridspec

import pandas as pd  ###To treat the csv files

from matplotlib import ticker

formatter = ticker.ScalarFormatter(useMathText=True)
formatter.set_scientific(True)
formatter.set_powerlimits((-3, 3))

###Defined Colors optimized for color-blind people:
Orange = [230 / 255, 159 / 255, 0 / 255]
SkyBlue = [86 / 255, 180 / 255, 233 / 255]
BluishGreen = [0 / 255, 158 / 255, 115 / 255]
Yellow = [240 / 255, 228 / 255, 66 / 255]
Blue = [0 / 255, 114 / 255, 178 / 255]
Vermillion = [213 / 255, 94 / 255, 0 / 255]
ReddishPurple = [204 / 255, 121 / 255, 167 / 255]
####//////////////////////////////////////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\####
####///~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\\\####
####///~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ MAIN PART OF THE CODE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\\\####
####///~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\\\####
####//////////////////////////////////////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\####
############################################################
# Define function to check whether node is above cavity ####
############################################################
def IsAboveCavity(x,y):
    xc = 948003
    Rx = 22
    yc = 2105064
    Ry = 54
    metric = ((xc - x)/Rx)**2 + ((yc - y)/Ry)**2
    out = metric <= 1 ###LOGICAL: True if above cavity, False otherwise
    return out

################################################################################
# prepare plots #####
################################################################################
# global figure pimping
# plt.rc('text', usetex=True)
plt.rc('xtick', labelsize=22)
plt.rc('ytick', labelsize=22)
plt.rc('axes', labelsize=22)
plt.rc('legend', fontsize=26)

nrows = 4
ncols = 2
fig, axes = plt.subplots(nrows, ncols, figsize=(80, 40), sharex= True , sharey= True)
# Common ylabel for all subplots
fig.text(0.005, 0.5, r'Mean $V_z$ [mm/d]', va='center', rotation='vertical', fontsize=24)
# pimp style
for i in range(nrows): ###
    for j in range(ncols):
        ax = axes[i,j]
        ax.tick_params(labelsize=18)  # fontsize of the tick labels
        ax.set_ylim([-11.0, 2.0])
        # ax.grid(True)
        ax.yaxis.set_major_formatter(formatter)
        ax.xaxis.set_major_locator(MonthLocator(interval=3))
        ax.xaxis.set_minor_locator(MonthLocator(interval=1))
        ax.xaxis.set_major_formatter(DateFormatter('%Y-%m'))
        fig.autofmt_xdate()
        ###plot one horizontal line for Vz=0
        ax.axhline(y=0.0, color= 'k', linestyle = ':', linewidth=1)
        ###Shade periods corresponding to pumping (26/08/2010 to 15/10/10; 28/09/11 to14/10/11; 23/09/12 to 09/10/12)
        ax.axvspan(date(int(2010), 8, 26) ,date(int(2010), 10, 15) , alpha=0.3, color='grey')
        ax.axvspan(date(int(2011), 9, 28) ,date(int(2011), 10, 14) , alpha=0.3, color='grey')
        ax.axvspan(date(int(2012), 9, 23) ,date(int(2012), 10, 9) , alpha=0.3, color='grey')



################################################################################
# Make the plots #####
################################################################################
if __name__ == "__main__":
    ####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####
    ####~~~~~~~~~~~   LOAD ALL DATA AND PLOT DIRECTLY IN THIS PART OF THE CODE   ~~~~~~~~~~~~~####
    ####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####
    ### General parameter of simu
    RefStartDate = date(int(2011), 7, 1) ##Every days are numbered relative to 1st July 2011 (day 1)
    StartDayNumber_Pumping2010 = "-315"
    StartYear_Pumping2010 = "2010"
    StartDate_Pumping2010 = RefStartDate + timedelta(days=int(StartDayNumber_Pumping2010))

    #########################################
    #### SOME PARAMETERS FOR THE PLOTS:  ####
    #########################################
    ###Load file containing all tested parameter sets for damage model
    Pathroot_ParamSet = Path('/home/brondexj/BETTIK/TeteRousse/MyTeterousse_GeoGag/.')
    Filename_ParamSet = 'Dam_Sigmath_B_Lambdah.IN'
    Col_Names_ParamSet = ['Sigmath', 'B', 'Lambdah']
    Data_ParamSet = pd.read_csv(Pathroot_ParamSet.joinpath(Filename_ParamSet), names=Col_Names_ParamSet, delim_whitespace=True)
    ###Add colums to Data_ParamSet corresponding to parameters converted in file names (convert to string and remove dot)
    Data_ParamSet['Sigmath_name']=Data_ParamSet['Sigmath'].astype(str).replace('\.', '',regex=True)
    Data_ParamSet['B_name']=Data_ParamSet['B'].astype(str).replace('\.', '',regex=True)
    Data_ParamSet['Lambdah_name']=Data_ParamSet['Lambdah'].astype(str).replace('\.', '',regex=True)
    ###List of tested lambdah values
    Lambdah_list = [0.0]
    ###Corresponding linestyle
    LineStyle_List = ['-']
    ###Alternative choice for colors: a discrete colormap that is colorblind-friendly
    ###In this order: blue/orange/green/pink/brown/purple/grey/red/yellow
    ColorBlindFriendly_List = ['#377eb8', '#ff7f00', '#4daf4a','#f781bf', '#a65628', '#984ea3','#999999', '#e41a1c', '#dede00']
    ### Load files containing surface output of the simulations:
    Pathroot_SimuOutput = Path('/home/brondexj/BETTIK/TeteRousse/MyTeterousse_GeoGag/ScalarOutput/.')
    ### Step in the full process from 2010 to 2014
    Step_List = ['P2010', 'Refill20102011', 'P2011', 'Refill20112012', 'P2012', 'Refill20122013']
    StepTsp_List =[1, 5, 1, 5, 1, 5] ##Timestep size (in days) corresponding to simulation step (Step_List)
    StepOut_List =[5, 30, 5, 30, 5, 30] ##Output interval (in days) corresponding to simulation step (Step_List)
    ### List of tested damage criterion and corresponding name
    Crit_List = ['MPS', 'Ha','vM', 'Co']
    Crit_Names = ['Max Principal Stress', 'Hayurst','von Mises', 'Coulomb']
    ### Where pressure applies: cavity only, restricted to a conduit, everywhere above cold/temperate transition
    case = 'PCav' #, 'PRest', 'PNotRest'
    ###NAMES OF OUTPUT DATA (same order as in dat.names file)
    Col_Names = ['Timestep', 'DayOfSimu', 'BC', 'NodeNumber', 'X', 'Y', 'Z', 'U', 'V', 'W', 'SigmaI', 'Damage', 'Chi']
    Col_Names_NoD_Tmap = ['Timestep', 'DayOfSimu', 'BC', 'NodeNumber', 'X', 'Y', 'Z', 'U', 'V', 'W', 'Pressure', 'SigmaIII','SigmaII', 'SigmaI', 'Temperature', 'Sxx', 'Syy', 'Szz', 'Sxy', 'Sxz', 'Syz' ]  ##For the no damage (ref) case
    Col_Names_NoD_Tfix = ['Timestep', 'DayOfSimu', 'BC', 'NodeNumber', 'X', 'Y', 'Z', 'U', 'V', 'W', 'Pressure', 'SigmaIII','SigmaII', 'SigmaI']

    ###################################
    #### START TREATMENT OF DATA:  ####
    ###################################
    ###FIRST OF ALL: deal with the No Damage case that will be used as reference for all damage combinations
    ###We create a single dataframe for all steps
    Data_Simu_NoD = pd.DataFrame()  ##Initialize an empty dataframe
    for i, (Step, StepTsp, StepOut) in enumerate(zip(Step_List, StepTsp_List, StepOut_List)):
        filename = 'SurfaceOutput_Prono_NoD_{}_{}_tsp{}d_Out{}d_Tmap_.dat'.format(case, Step, str(StepTsp), str(StepOut))
        print('Opening file:', filename)
        Data_Simu_NoD_tmp = pd.read_csv(Pathroot_SimuOutput.joinpath(filename), names=Col_Names_NoD_Tmap, delim_whitespace=True)
        ###Drop duplicate lines (bug of SaveLine solver)
        Data_Simu_NoD_tmp.drop_duplicates(inplace=True)
        ###Set Day of Simu counting from begining of simu corresponding to step Pump2010 (Step 1)
        if i == 0:
            Data_Simu_NoD_tmp['DayOfSimu'] = Data_Simu_NoD_tmp['DayOfSimu'] * StepTsp
        else:
            Data_Simu_NoD_tmp['DayOfSimu'] = np.max(Data_Simu_NoD['DayOfSimu']) + Data_Simu_NoD_tmp['DayOfSimu'] * StepTsp
        data = [Data_Simu_NoD, Data_Simu_NoD_tmp]
        Data_Simu_NoD = pd.concat(data, ignore_index=True)
    ###Add a column of boolean to dataset depending on wether node is above cavity or not
    Data_Simu_NoD['IsAboveCavity'] = IsAboveCavity(Data_Simu_NoD['X'], Data_Simu_NoD['Y'])
    ###Store Node above cavity in a separated dataframe
    Data_Simu_NoD_AboveCavity = Data_Simu_NoD[Data_Simu_NoD['IsAboveCavity']]
    ###For each timestep calculate mean, min, max of vertical velocity above cavity
    Date = []
    MeanW_NoD = []
    MinW_NoD = []
    MaxW_NoD = []
    ### We want to produce a list corresponding to all simulation days at which we have an output
    SimuDays = Data_Simu_NoD['DayOfSimu'].drop_duplicates().reset_index(drop=True)
    for k in range(len(SimuDays)):
        day = int(SimuDays.iloc[k])
        ###Convert timestep in terms of date
        Date.append(StartDate_Pumping2010 + timedelta(days=day - 1))
        ###Calculate mean/min/max in mm/days
        MeanW_NoD.append(np.mean(Data_Simu_NoD_AboveCavity[Data_Simu_NoD_AboveCavity['DayOfSimu'] == day]['W']) * (1000 / 365.25))
        MinW_NoD.append(np.min(Data_Simu_NoD_AboveCavity[Data_Simu_NoD_AboveCavity['DayOfSimu'] == day]['W']) * (1000 / 365.25))
        MaxW_NoD.append(np.max(Data_Simu_NoD_AboveCavity[Data_Simu_NoD_AboveCavity['DayOfSimu'] == day]['W']) * (1000 / 365.25))
    ###~~~~~~~~~~~~~SAME FOR SIMULATION CONSIDERING TEMPERATE ICE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ###FIRST OF ALL: deal with the No Damage case that will be used as reference for all damage combinations
    ###We create a single dataframe for all steps
    Data_Simu_NoD_Temperate = pd.DataFrame()  ##Initialize an empty dataframe
    for i, (Step, StepTsp, StepOut) in enumerate(zip(Step_List, StepTsp_List, StepOut_List)):
        filename = 'SurfaceOutput_Prono_NoD_{}_{}_tsp{}d_Out{}d_Temperate_.dat'.format(case, Step, str(StepTsp), str(StepOut))
        print('Opening file:', filename)
        Data_Simu_NoD_tmp_Temperate = pd.read_csv(Pathroot_SimuOutput.joinpath(filename), names=Col_Names_NoD_Tfix,delim_whitespace=True)
        ###Drop duplicate lines (bug of SaveLine solver)
        Data_Simu_NoD_tmp_Temperate.drop_duplicates(inplace=True)
        ###Set Day of Simu counting from begining of simu corresponding to step Pump2010 (Step 1)
        if i == 0:
            Data_Simu_NoD_tmp_Temperate['DayOfSimu'] = Data_Simu_NoD_tmp_Temperate['DayOfSimu'] * StepTsp
        else:
            Data_Simu_NoD_tmp_Temperate['DayOfSimu'] = np.max(Data_Simu_NoD_Temperate['DayOfSimu']) + Data_Simu_NoD_tmp_Temperate['DayOfSimu'] * StepTsp
        data_Temperate = [Data_Simu_NoD_Temperate, Data_Simu_NoD_tmp_Temperate]
        Data_Simu_NoD_Temperate = pd.concat(data_Temperate, ignore_index=True)
    ###Add a column of boolean to dataset depending on wether node is above cavity or not
    Data_Simu_NoD_Temperate['IsAboveCavity'] = IsAboveCavity(Data_Simu_NoD_Temperate['X'], Data_Simu_NoD_Temperate['Y'])
    ###Store Node above cavity in a separated dataframe
    Data_Simu_NoD_AboveCavity_Temperate = Data_Simu_NoD_Temperate[Data_Simu_NoD_Temperate['IsAboveCavity']]
    ###For each timestep calculate mean, min, max of vertical velocity above cavity
   # Date = []
    MeanW_NoD_Temperate = []
    MinW_NoD_Temperate = []
    MaxW_NoD_Temperate = []
    ### We want to produce a list corresponding to all simulation days at which we have an output
#    SimuDays = Data_Simu_NoD_Temperate['DayOfSimu'].drop_duplicates().reset_index(drop=True)
    for k in range(len(SimuDays)):
        day = int(SimuDays.iloc[k])
        ###Convert timestep in terms of date
       # Date.append(StartDate_Pumping2010 + timedelta(days=day - 1))
        ###Calculate mean/min/max in mm/days
        MeanW_NoD_Temperate.append(
            np.mean(Data_Simu_NoD_AboveCavity_Temperate[Data_Simu_NoD_AboveCavity_Temperate['DayOfSimu'] == day]['W']) * (1000 / 365.25))
        MinW_NoD_Temperate.append(
            np.min(Data_Simu_NoD_AboveCavity_Temperate[Data_Simu_NoD_AboveCavity_Temperate['DayOfSimu'] == day]['W']) * (1000 / 365.25))
        MaxW_NoD_Temperate.append(
            np.max(Data_Simu_NoD_AboveCavity_Temperate[Data_Simu_NoD_AboveCavity_Temperate['DayOfSimu'] == day]['W']) * (1000 / 365.25))
    ###~~~~~~~~~~~~~SAME FOR SIMULATION CONSIDERING ICE AT T=-2°C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ###We create a single dataframe for all steps
    Data_Simu_NoD_Tminus2 = pd.DataFrame()  ##Initialize an empty dataframe
    for i, (Step, StepTsp, StepOut) in enumerate(zip(Step_List, StepTsp_List, StepOut_List)):
        filename = 'SurfaceOutput_Prono_NoD_{}_{}_tsp{}d_Out{}d_Tminus2_.dat'.format(case, Step, str(StepTsp),
                                                                                       str(StepOut))
        print('Opening file:', filename)
        Data_Simu_NoD_tmp_Tminus2 = pd.read_csv(Pathroot_SimuOutput.joinpath(filename), names=Col_Names_NoD_Tfix,
                                                  delim_whitespace=True)
        ###Drop duplicate lines (bug of SaveLine solver)
        Data_Simu_NoD_tmp_Tminus2.drop_duplicates(inplace=True)
        ###Set Day of Simu counting from begining of simu corresponding to step Pump2010 (Step 1)
        if i == 0:
            Data_Simu_NoD_tmp_Tminus2['DayOfSimu'] = Data_Simu_NoD_tmp_Tminus2['DayOfSimu'] * StepTsp
        else:
            Data_Simu_NoD_tmp_Tminus2['DayOfSimu'] = np.max(Data_Simu_NoD_Tminus2['DayOfSimu']) + \
                                                       Data_Simu_NoD_tmp_Tminus2['DayOfSimu'] * StepTsp
        data_Tminus2 = [Data_Simu_NoD_Tminus2, Data_Simu_NoD_tmp_Tminus2]
        Data_Simu_NoD_Tminus2 = pd.concat(data_Tminus2, ignore_index=True)
    ###Add a column of boolean to dataset depending on wether node is above cavity or not
    Data_Simu_NoD_Tminus2['IsAboveCavity'] = IsAboveCavity(Data_Simu_NoD_Tminus2['X'],
                                                             Data_Simu_NoD_Tminus2['Y'])
    ###Store Node above cavity in a separated dataframe
    Data_Simu_NoD_AboveCavity_Tminus2 = Data_Simu_NoD_Tminus2[Data_Simu_NoD_Tminus2['IsAboveCavity']]
    ###For each timestep calculate mean, min, max of vertical velocity above cavity
    #Date = []
    MeanW_NoD_Tminus2 = []
    MinW_NoD_Tminus2 = []
    MaxW_NoD_Tminus2 = []
    ### We want to produce a list corresponding to all simulation days at which we have an output
   # SimuDays = Data_Simu_NoD_Tminus2['DayOfSimu'].drop_duplicates().reset_index(drop=True)
    for k in range(len(SimuDays)):
        day = int(SimuDays.iloc[k])
        ###Convert timestep in terms of date
      #  Date.append(StartDate_Pumping2010 + timedelta(days=day - 1))
        ###Calculate mean/min/max in mm/days
        MeanW_NoD_Tminus2.append(
            np.mean(Data_Simu_NoD_AboveCavity_Tminus2[Data_Simu_NoD_AboveCavity_Tminus2['DayOfSimu'] == day][
                        'W']) * (1000 / 365.25))
        MinW_NoD_Tminus2.append(
            np.min(Data_Simu_NoD_AboveCavity_Tminus2[Data_Simu_NoD_AboveCavity_Tminus2['DayOfSimu'] == day][
                       'W']) * (1000 / 365.25))
        MaxW_NoD_Tminus2.append(
            np.max(Data_Simu_NoD_AboveCavity_Tminus2[Data_Simu_NoD_AboveCavity_Tminus2['DayOfSimu'] == day][
                       'W']) * (1000 / 365.25))

    # # ### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
    # # ### NOW START LOOP over each damage parameters combination (one subplot per combination)
    for k in range((Data_ParamSet['Lambdah']==0.0).sum()):
        ### Get the corresponding subplot
        i = k % nrows
        j = k // nrows
        ax = axes[i,j]
        # ax.set_title('Combi {}'.format(k+1), fontsize=20, weight='bold') ###Easier to add afterwards on inkscape
        ###NOW DO THE LOOP ON EACH DAMAGE CRITERION
        for l,(Crit,Crit_Name) in enumerate(zip(Crit_List,Crit_Names)):##range(len(Data_ParamSet)):
            print('Considering Combi',k+1,'and criterion ',Crit_Name)
            ###We create a single dataframe for all steps
            Data_Simu_D = pd.DataFrame()  ##Initialize an empty dataframe
            for i, (Step, StepTsp, StepOut) in enumerate(zip(Step_List, StepTsp_List, StepOut_List)):
                filename = 'SurfaceOutput_Prono_B{}_Sth{}_Lh{}_{}_{}_tsp{}d_Out{}d_Crit{}_Tmap_.dat'.format(Data_ParamSet['B_name'][k], Data_ParamSet['Sigmath_name'][k], Data_ParamSet['Lambdah_name'][k], case, Step, str(StepTsp), str(StepOut), Crit)
                print('Opening file:', filename)
                Data_Simu_D_tmp = pd.read_csv(Pathroot_SimuOutput.joinpath(filename), names=Col_Names, delim_whitespace=True)
                ###Drop duplicate lines (bug of SaveLine solver
                Data_Simu_D_tmp.drop_duplicates(inplace=True)
                ###Set Day of Simu counting from begining of simu corresponding to step Pump2010 (Step 1)
                if i == 0:
                    Data_Simu_D_tmp['DayOfSimu'] = Data_Simu_D_tmp['DayOfSimu'] * StepTsp
                else:
                    Data_Simu_D_tmp['DayOfSimu'] = np.max(Data_Simu_D['DayOfSimu']) + Data_Simu_D_tmp['DayOfSimu'] * StepTsp
                data_1 = [Data_Simu_D, Data_Simu_D_tmp]
                Data_Simu_D = pd.concat(data_1, ignore_index=True)
            ###Add a column of boolean to dataset depending on wether node is above cavity or not
            Data_Simu_D['IsAboveCavity']=IsAboveCavity(Data_Simu_D['X'],Data_Simu_D['Y'])
            ###Store Node above cavity in a separated dataframe
            Data_Simu_D_AboveCavity=Data_Simu_D[Data_Simu_D['IsAboveCavity']]
            ###For each timestep calculate mean, min, max of vertical velocity above cavity
            Date = []
            MeanW_D = []
            MinW_D = []
            MaxW_D = []
            ### We want to product a list corresponding to all simulation days at which we have an output
            SimuDays = Data_Simu_D['DayOfSimu'].drop_duplicates().reset_index(drop=True)
            for m in range(len(SimuDays)):
                day = int(SimuDays.iloc[m])
                ###Convert timestep in terms of date
                Date.append(StartDate_Pumping2010 + timedelta(days=day-1))
                ###Calculate mean/min/max in mm/days
                MeanW_D.append(np.mean(Data_Simu_D_AboveCavity[Data_Simu_D_AboveCavity['DayOfSimu']==day]['W']) * (1000/365.25))
                MinW_D.append(np.min(Data_Simu_D_AboveCavity[Data_Simu_D_AboveCavity['DayOfSimu']==day]['W']) * (1000/365.25))
                MaxW_D.append(np.max(Data_Simu_D_AboveCavity[Data_Simu_D_AboveCavity['DayOfSimu']==day]['W']) * (1000/365.25))
            ##Plot MeanW for the case with damage on corresponding subplot
            Color = ColorBlindFriendly_List[l]
            ax.plot_date(Date, MeanW_D, color= Color, linestyle = '-', linewidth=1.8, marker='None', xdate=True)
        ###Plot MeanW for the case with no damage with temperate ice on every subplot
        ax.plot_date(Date, MeanW_NoD_Temperate, color= 'k', linestyle = ':', linewidth=1.7, marker='None', xdate=True)
        ###Plot MeanW for the case with no damage with ice at T=-2°C on every subplot
        ax.plot_date(Date, MeanW_NoD_Tminus2, color= 'k', linestyle = '--', linewidth=1.7, marker='None', xdate=True)
        ###Plot MeanW for the case with no damage (ref case) on every subplot
        ax.plot_date(Date, MeanW_NoD, color= 'k', linestyle = '-', linewidth=2, marker='None', xdate=True)

    #### DUMMY plot for legend
    ax=axes[0,0]
    ax.plot(np.nan, np.nan, label=r'No Damage, T = f(x,z)', color='k', linewidth=3, linestyle='-')
    ax.plot(np.nan, np.nan, label=r'No Damage, T = 0°C', color='k', linewidth=3, linestyle=':')
    ax.plot(np.nan, np.nan, label=r'No Damage, T = -2°C', color='k', linewidth=3, linestyle='--')
    for i,(Crit,Crit_Name) in enumerate(zip(Crit_List,Crit_Names)):
        ax.plot(np.nan, np.nan, label=r'{}'.format(Crit_Name), color=ColorBlindFriendly_List[i], linewidth=3, linestyle='-')
    fig.legend(loc='lower center', bbox_to_anchor=(0.5, 0.005), fancybox=True, shadow=True, fontsize=16, ncol=4)
    # Adjust layout to make room for the legend
    plt.tight_layout(rect=[0.05, 0.17, 0.96, 0.96],h_pad=6.0, w_pad=5.0)

    #########################################
    ####   NOW OBSERVATIONS AT STAKES    ####
    #########################################
    ####Load data corresponding to observations at stake
    Pathroot_ObsData = Path('/home/brondexj/BETTIK/TeteRousse/MyTeterousse_GeoGag/PostProcessing/Raw_Data_VeloAtStakes_CORRECTED/.')
    ###List of all files containing observationnal data
    ObsDataName_List = ['coord14092010','coord19092010', 'coord20092010', 'coord21092010', 'coord22092010',
                        'coord23092010', 'coord29092010','coord06102010',
                        'coord25052011', 'coord15062011', 'coord01072011', 'coord15072011', 'coord10082011',
                        'coord31082011', 'coord09092011', 'coord23092011', 'coord28092011', 'coord06102011',
                        'coord21102011', 'coord17112011',
                        'coord30082011-pieubois', 'coord09092011-pieubois', 'coord28092011-pieubois',
                        'coord13072012', 'coord25072012', 'coord09082012', 'coord22082012', 'coord05092012',
                        'coord20092012', 'coord28092012', 'coord01102012', 'coord03102012']
    ###For each year, list of stakes that are not above of the cavity (i.e. out of the AboveCavity mask)
    StakeNotAboveCavity_2010 = [1, 11, 12, 17, 18, 22, 27] ###Stakes 10, 13 are also out but close
    StakeNotAboveCavity_2011 = [4, 5, 6, 15, 16, 20, 21, 24, 25, 29, 30] ###Stakes 7, 11, 26, 28 also out but close
    StakeNotAboveCavity_2011_pieubois = [81, 82, 83, 87, 92, 93, 99] ###Stakes 80 also out but close
    StakeNotAboveCavity_2012 = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 16, 17, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44] ##Stake 18 also out but close
    ###NAMES COLUMS of obs data (same order as in dat.names file)
    Col_Names_ObsData = ['Stake', 'X', 'Y', 'Z']  ##obs data file contains stake number, coords x, y and z
    ### PRESCRIBE THE DISPLACEMENTS THAT WE WANT TO CALCULATE FROM AVAILABLE OBS DATA
    Displacement_List = ['disp14092010-19092010', 'disp19092010-23092010', 'disp23092010-29092010', 'disp29092010-06102010',
                         'disp25052011-15062011', 'disp15062011-01072011', 'disp01072011-15072011', 'disp15072011-10082011',
                         'disp10082011-31082011', 'disp31082011-09092011', 'disp09092011-23092011', 'disp23092011-06102011',
                         'disp06102011-21102011', 'disp21102011-17112011',
    #                     'disp30082011-09092011-pieubois', 'disp09092011-28092011-pieubois',
                         'disp13072012-25072012', 'disp25072012-05092012', 'disp05092012-20092012', 'disp20092012-28092012',
                         'disp28092012-03102012'
                         ]
    ###DO THE LOOP ON ALL OBS DATA FILES TO BUILD A DATAFRAME STORING ALL OBS
    ###We build a single dataframe containing all coord for all days
    Coords_df = pd.DataFrame()
    for i, ObsDataName in enumerate(ObsDataName_List):
        date_part = ObsDataName[5:13]
        day = int(date_part[:2])
        month = int(date_part[2:4])
        year = int(date_part[4:])
        pieubois = 'pieubois' in ObsDataName
        Path_To_ObsFile = 'topo{}/{}.dat'.format(str(year), ObsDataName)
        ObsData_tmp = pd.read_csv(Pathroot_ObsData.joinpath(Path_To_ObsFile), names=Col_Names_ObsData, delim_whitespace=True)
        ##Remove stakes that are not above cavity
        for l in range(len(ObsData_tmp)):
            StakeNo = ObsData_tmp['Stake'][l]
            if ((year == 2010) and (StakeNo in StakeNotAboveCavity_2010)):
                ObsData_tmp.drop([l], inplace=True)
            elif ((year == 2011) and (not pieubois) and (StakeNo in StakeNotAboveCavity_2011)):
                ObsData_tmp.drop([l], inplace=True)
            elif ((year == 2011) and (pieubois) and (StakeNo in StakeNotAboveCavity_2011_pieubois)):
                ObsData_tmp.drop([l], inplace=True)
            elif ((year == 2012) and (StakeNo in StakeNotAboveCavity_2012)):
                ObsData_tmp.drop([l], inplace=True)
        Coords_df_tmp = ObsData_tmp.copy()
        Coords_df_tmp['day']=day
        Coords_df_tmp['month']=month
        Coords_df_tmp['year']=year
        Coords_df_tmp['pieubois']=pieubois
        data = [Coords_df, Coords_df_tmp]
        Coords_df = pd.concat(data, ignore_index=True)

    ###NOW DO THE LOOP ON ALL PRESCRIBED DISPLACEMENTS AND DEDUCE VELOCITY
    ###We create lists for Date, Wmean, Wmin and Wmax
    Obs_Date = []
    Obs_MeanW = []
    Obs_MinW = []
    Obs_MaxW = []
    for i, Disp_Name in enumerate(Displacement_List):
        pieubois_disp = 'pieubois' in Disp_Name

        startdate_part = Disp_Name[4:12]
        startday = int(startdate_part[:2])
        startmonth = int(startdate_part[2:4])
        startyear = int(startdate_part[4:])
        startdate = datetime(startyear, startmonth, startday, 0, 0)

        finishdate_part = Disp_Name[13:21]
        finishday = int(finishdate_part[:2])
        finishmonth = int(finishdate_part[2:4])
        finishyear = int(finishdate_part[4:])
        finishdate = datetime(finishyear, finishmonth, finishday, 0, 0)

        NumberOfDays = (finishdate - startdate).days
        MeanDay = startdate + (finishdate - startdate)/2
        Obs_Date.append(MeanDay)
        ###Get altitude of stakes for start and finish dates
        StartDate_Coord_df = Coords_df[(Coords_df['day']==startday) & (Coords_df['month']==startmonth) & (Coords_df['year']==startyear) & (Coords_df['pieubois']==pieubois_disp)]
        FinishDate_Coord_df = Coords_df[(Coords_df['day']==finishday) & (Coords_df['month']==finishmonth) & (Coords_df['year']==finishyear) & (Coords_df['pieubois']==pieubois_disp)]
        ###Calculate displacements between two dates on a stake by stake basis
        Disp_df = pd.DataFrame()
        # Filter for stakes present in both DataFrames
        common_stakes = StartDate_Coord_df[StartDate_Coord_df["Stake"].isin(FinishDate_Coord_df["Stake"])]
        # Calculate displacements
        for _, row in common_stakes.iterrows():
            stake = row['Stake']
            start_z = row['Z']
            finish_z = FinishDate_Coord_df[FinishDate_Coord_df['Stake'] == stake]['Z'].values[0]
            # Append to Disp_df
            Disp_df = pd.concat([Disp_df, pd.DataFrame({'Stake': [stake], 'Uz': [finish_z - start_z]}),], ignore_index=True,)
        ###Calculate and store mean, min and max W over the period
        print('Displacement period:', Disp_Name, 'Retained displacement date:', MeanDay)
        MeanW = np.mean(Disp_df['Uz'])*1000 / NumberOfDays
        Obs_MeanW.append(MeanW)
        MinW = np.min(Disp_df['Uz'])*1000 / NumberOfDays
        Obs_MinW.append(MinW)
        MaxW = np.max(Disp_df['Uz'])*1000 / NumberOfDays
        Obs_MaxW.append(MaxW)

    ###Now plot point (mean) with error bar min/max on each subplot
    for ax in axes.reshape(-1):
        ax.errorbar(Obs_Date, Obs_MeanW,yerr=[abs(np.subtract(Obs_MeanW, Obs_MinW)), abs(np.subtract(Obs_MeanW, Obs_MaxW))], fmt='o',
                    markerfacecolor='w', capsize=3)

    plt.show()
