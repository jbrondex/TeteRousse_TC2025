"""
@author: jbrondex

Description:
------------
This file aims at plotting the displacement field above the cavity computed between two prescribed dates to compare displacement predicted by
elastic and viscous constitutive laws to what has been actually observed

"""

################################################################################
# import #####
################################################################################
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm, ListedColormap
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
from scipy.interpolate import griddata

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
##############################
# Define useful functions ####
##############################
def IsAboveCavity(x,y): ####Define function to check whether node is above cavity
    xc = 948003
    Rx = 22
    yc = 2105064
    Ry = 54
    metric = ((xc - x)/Rx)**2 + ((yc - y)/Ry)**2
    out = metric <= 1 ###LOGICAL: True if above cavity, False otherwise
    return out

def Interpolate_field(df, field_name, x, y): ###Returned field interpolated from mesh grid to point (x,y) of transect
    # x, y, and z coordinates of the dataframe
    xi = df['X'].values
    yi = df['Y'].values
    field = df[field_name].values
    # Interpolate the value of field
    result = griddata((xi, yi), field, (x, y), rescale=True)
    return result

##Function below is used to sort point of contour clockwise
def sort_clockwise(list_of_xy_coords):
    cx, cy = list_of_xy_coords.mean(0)
    x, y = list_of_xy_coords.T
    angles = np.arctan2(x-cx, y-cy)
    indices = np.argsort(angles)
    return list_of_xy_coords[indices]


# Interpolation function for day of observations that are not among day of simu
def interpolate_at_date(df, day_col, obs_day):
    # Ensure the DataFrame is sorted by the day column
    df = df.sort_values(by=day_col).reset_index(drop=True)

    # Find the two closest neighbors (days right before and after the target_day)
    before = df[df[day_col] < obs_day].iloc[-1]
    after = df[df[day_col] > obs_day].iloc[0]

    # Linear interpolation for each node column
    interpolated_values = {}
    for col in df.columns:
        if col == day_col:  # Skip the day column
            continue
        # Linear interpolation formula
        interpolated_values[col] = before[col] + (after[col] - before[col]) * (
                (obs_day - before[day_col]) / (after[day_col] - before[day_col])
        )

    # Combine interpolated values with the target day
    interpolated_row = {day_col: obs_day, **interpolated_values}
    return interpolated_row


################################################################################
# prepare plots #####
################################################################################
# global figure pimping
# plt.rc('text', usetex=True)
plt.rc('xtick', labelsize=22)
plt.rc('ytick', labelsize=22)
plt.rc('axes', labelsize=22)
plt.rc('legend', fontsize=26)

nrows = 3 ###Ux, Uy, Uz
ncols = 3 ###elastic, viscous, observed
fig, axes = plt.subplots(nrows, ncols, figsize=(15, 12), sharex= True , sharey= True, constrained_layout=True)
# Common xlabel for all subplots
fig.text(0.5, 0.005, r'X [km]', va='center', rotation='horizontal', fontsize=24)
# Common ylabel for all subplots
fig.text(0.005, 0.5, r'Y [km]', va='center', rotation='vertical', fontsize=24)
# pimp style
for i in range(nrows): ###
    for j in range(ncols):
        ax = axes[i,j]
        ax.tick_params(labelsize=18)  # fontsize of the tick labels
        ax.grid(True)
        ax.grid(alpha=0.5)
        ###Force x and y axis to same scale
        ax.set_aspect('equal', adjustable='box')


################################################################################
# Make the plots #####
################################################################################
if __name__ == "__main__":
    ####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####
    ####~~~~~~~~~~~   LOAD ALL DATA AND PLOT DIRECTLY IN THIS PART OF THE CODE   ~~~~~~~~~~~~~####
    ####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####
    ### General parameter regarding dates of interest
    RefStartDate = date(int(2011), 7, 1) ##Every days are numbered relative to 1st July 2011 (day 1)
    StartDate_str = '14092010'  ###start date for displacement computation
    FinishDate_str = '06102010' ###finish date for displacement computation
    Sequence = 'Pumping2010' ###To which sequence corresponds the chosen dates ? Required to set proper legend parameters
    # StartDate_str = '09092011'  ###start date for displacement computation
    # FinishDate_str = '28092011' ###finish date for displacement computation
    # Sequence = 'Refill20102011'
    # StartDate_str = '28092011'  ###start date for displacement computation
    # FinishDate_str = '21102011' ###finish date for displacement computation
    # Sequence = 'Pumping2011'
    ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####
    ### Parameter for the plot
    # shading
    cmap = 'coolwarm'
    color_resolution = 4  ##number of color shades per level ticks
    norm = None ##by default (to keep only positive or negative side of colorbar)
    if Sequence == 'Pumping2010':
        dx_min_el = -0.16
        dx_max_el = 0.161
        lev_ticks_range_dx_el = 0.08 ##one tick every ?
        norm_dx_el = norm
        dx_min_obs = -14.0
        dx_max_obs = 14.1
        lev_ticks_range_dx_obs = 7 ##one tick every ?
        norm_dx_obs = norm
        #####
        dy_min_el = -0.10
        dy_max_el = 0.11
        lev_ticks_range_dy_el = 0.05 ##one tick every ?
        norm_dy_el = norm
        dy_min_obs =-10.0
        dy_max_obs =10.1
        lev_ticks_range_dy_obs = 5 ##one tick every ?
        norm_dy_obs = norm
        #####
        dz_min_el =-0.24
        dz_max_el =0.001
        lev_ticks_range_dz_el = 0.08 ##one tick every ?
        norm_dz_el = TwoSlopeNorm(vmin=dz_min_el, vcenter=0.0, vmax=dz_max_el)
        dz_min_obs =-24.0
        dz_max_obs =0.1
        lev_ticks_range_dz_obs = 8 ##one tick every ?
        norm_dz_obs = TwoSlopeNorm(vmin=dz_min_obs, vcenter=0.0, vmax=dz_max_obs)
    ################################################################
    elif Sequence == 'Refill20102011':
        dx_min_el = -0.05
        dx_max_el = 0.051
        lev_ticks_range_dx_el = 0.025 ##one tick every ?
        norm_dx_el = norm
        dx_min_obs = -6.0
        dx_max_obs = 0.1
        lev_ticks_range_dx_obs = 1.5 ##one tick every ?
        norm_dx_obs = TwoSlopeNorm(vmin=dx_min_obs, vcenter=0.0, vmax=dx_max_obs)
        #####
        dy_min_el = -0.04
        dy_max_el = 0.041
        lev_ticks_range_dy_el = 0.02 ##one tick every ?
        norm_dy_el = norm
        dy_min_obs =-2.0
        dy_max_obs =2.1
        lev_ticks_range_dy_obs = 1 ##one tick every ?
        norm_dy_obs = norm
        #####
        dz_min_el =0.000
        dz_max_el =0.11
        lev_ticks_range_dz_el = 0.025*2 ##one tick every ?
        norm_dz_el = TwoSlopeNorm(vmin=dz_min_el-0.00001, vcenter=0.0, vmax=dz_max_el)
        dz_min_obs =-5
        dz_max_obs =5.1
        lev_ticks_range_dz_obs = 2.5 ##one tick every ?
        norm_dz_obs = norm
    ##############################################################################
    elif Sequence == 'Pumping2011':
        dx_min_el = -0.30
        dx_max_el = 0.31
        lev_ticks_range_dx_el = 0.15 ##one tick every ?
        norm_dx_el = norm
        dx_min_obs = -10.0
        dx_max_obs = 10.1
        lev_ticks_range_dx_obs = 5 ##one tick every ?
        norm_dx_obs = norm
        #####
        dy_min_el = -0.25
        dy_max_el = 0.26
        lev_ticks_range_dy_el = 0.125 ##one tick every ?
        norm_dy_el = norm
        dy_min_obs =-8.0
        dy_max_obs =8.1
        lev_ticks_range_dy_obs = 4 ##one tick every ?
        norm_dy_obs = norm
        #####
        dz_min_el =-0.6
        dz_max_el =0.001
        lev_ticks_range_dz_el = 0.2 ##one tick every ?
        norm_dz_el = TwoSlopeNorm(vmin=dz_min_el, vcenter=0.0, vmax=dz_max_el)
        dz_min_obs =-20.0
        dz_max_obs =0.1
        lev_ticks_range_dz_obs = 5 ##one tick every ?
        norm_dz_obs = TwoSlopeNorm(vmin=dz_min_obs, vcenter=0.0, vmax=dz_max_obs)


    ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####
    ### Paths to files to open                ####
    ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####
    ###Path to folder where to find data of interest
    Pathroot = Path('/home/brondexj/BETTIK/TeteRousse/MyTeterousse_GeoGag/.')
    ###Name of columns of output files
    Col_Names_Viscous = ['Timestep', 'Timestep2', 'BC', 'NodeNumber', 'X', 'Y', 'Z', 'U', 'V', 'W', 'Pressure', 'SigmaIII','SigmaII', 'SigmaI', 'Temperature', 'Sxx', 'Syy', 'Szz', 'Sxy', 'Sxz', 'Syz' ]
    Col_Names_Elastic = ['Timestep', 'BC', 'NodeNumber', 'X', 'Y', 'Z', 'DispX', 'DispY', 'DispZ', 'VonMises', 'SigmaI','SigmaII', 'SigmaIII']
    Col_Names_ObsData = ['Stake', 'X', 'Y', 'Z']
    ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####
    ###Open Data corresponding to cavity contour from grounded mask (corrected)
    Pathroot_GM = Path('/home/brondexj/BETTIK/TeteRousse/MyTeterousse_GeoGag/ScalarOutput/.')
    Filename_GM = 'GroundedLine_CLEANED.dat'
    Col_Names_GM = ['Timestep', 'BC', 'NodeNumber', 'X', 'Y', 'Z', 'GM']
    ###load file as dataframe
    Df_GL = pd.read_csv(Pathroot_GM.joinpath(Filename_GM), names=Col_Names_GM, delim_whitespace=True)
    ###Sort coordinates of cavity contour clockwise
    xy_cavity = np.array([Df_GL['X'].values, Df_GL['Y'].values])
    xy_cavity = np.transpose(xy_cavity)
    xy_cavity_sorted = sort_clockwise(xy_cavity)
    ##close contour
    cavity_contour = np.vstack([xy_cavity_sorted, xy_cavity_sorted[0]])
    ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####
    ###Open Data corresponding to sonar contour (paper Gag 2011)
    Pathroot_GLsonar = Path('/home/brondexj/BETTIK/TeteRousse/MyTeterousse_GeoGag/Data/.')
    Filename_GLsonar = 'limite_cavity_ordered.xyz'
    Col_Names_GLsonar = [ 'X', 'Y', 'Z', 'dummy']
    ###load file as dataframe
    Df_GLsonar = pd.read_csv(Pathroot_GLsonar.joinpath(Filename_GLsonar), names=Col_Names_GLsonar, delim_whitespace=True)
    ##############################################
    ##############################################
    ##### LOAD DATA AND DO SOME PROCESSING #######
    ##############################################
    #################################
    ###load and process obs data ####
    #################################
    FileName_ObsFile_StartDay = './PostProcessing/Raw_Data_VeloAtStakes_CORRECTED/topo{}/coord{}.dat'.format(StartDate_str[4:], StartDate_str)
    FileName_ObsFile_FinishDay ='./PostProcessing/Raw_Data_VeloAtStakes_CORRECTED/topo{}/coord{}.dat'.format(FinishDate_str[4:], FinishDate_str)
    df_Obs_StartDay = pd.read_csv(Pathroot.joinpath(FileName_ObsFile_StartDay), names=Col_Names_ObsData, delim_whitespace=True)
    df_Obs_FinishDay = pd.read_csv(Pathroot.joinpath(FileName_ObsFile_FinishDay), names=Col_Names_ObsData, delim_whitespace=True)
    # Filter for stakes present in both DataFrames and merge them
    common_stakes = df_Obs_StartDay.merge(df_Obs_FinishDay,on="Stake",suffixes=("_Start", "_Finish"),)
    # Calculate displacements
    common_stakes["RelDispX"] = common_stakes["X_Finish"] - common_stakes["X_Start"]
    common_stakes["RelDispY"] = common_stakes["Y_Finish"] - common_stakes["Y_Start"]
    common_stakes["RelDispZ"] = common_stakes["Z_Finish"] - common_stakes["Z_Start"]
    # Create the new DataFrame with start day coordinates as reference
    df_Obs_DiffBetweenTwoDates = common_stakes[
        ["Stake", "X_Start", "Y_Start", "Z_Start", "RelDispX", "RelDispY", "RelDispZ"]].rename(
        columns={"X_Start": "X", "Y_Start": "Y", "Z_Start": "Z"}
    )
    # For some period, some stake seem to give very strange values (outliers ?): filter them out:
    if StartDate_str == '09092011' and FinishDate_str == '28092011': ##for this period stake 4 seems to be outlier
        df_Obs_DiffBetweenTwoDates = df_Obs_DiffBetweenTwoDates[~df_Obs_DiffBetweenTwoDates['Stake'].isin([4])]
    ###From this dataframe create a meshgrid on which every data will be interpolated for comparison
    ###Create refined regular grid and interpolate field over grid
    xmin = np.floor(np.min(df_Obs_DiffBetweenTwoDates['X'])) - 40
    xmax = xmin + 150
    ymin = np.floor(np.min(df_Obs_DiffBetweenTwoDates['Y'])) - 18
    ymax = ymin + 150
    x = np.arange(xmin, xmax, 1)
    y = np.arange(ymin, ymax, 1)
    X, Y = np.meshgrid(x, y)
    ###########################################################################
    ###load and process data for elastic simulations (two diag simulations)####
    ###########################################################################
    filename_StartDay = './ScalarOutput/SurfaceOutput_DiagElastic_Poisson03_E1GPa_CavityPressure_{}_NoDispLat_.dat'.format(StartDate_str)
    filename_FinishDay = './ScalarOutput/SurfaceOutput_DiagElastic_Poisson03_E1GPa_CavityPressure_{}_NoDispLat_.dat'.format(FinishDate_str)
    df_Elastic_StartDay = pd.read_csv(Pathroot.joinpath(filename_StartDay), names=Col_Names_Elastic, delim_whitespace=True)
    df_Elastic_FinishDay = pd.read_csv(Pathroot.joinpath(filename_FinishDay), names=Col_Names_Elastic, delim_whitespace=True)
    # Verify the first 3 columns are identical ['DayOfSimu', 'BC', 'NodeNumber'] (slight differences are seen for coordinates
    if (df_Elastic_StartDay.iloc[:, :3].equals(df_Elastic_FinishDay.iloc[:, :3])):
        print("The first 3 columns are identical.")
        # Keep the first 6 columns in the new DataFrame
        df_Elastic_DiffBetweenTwoDates = df_Elastic_StartDay.iloc[:, :3].copy()
        # Add columns for coordinates x, y and z
        df_Elastic_DiffBetweenTwoDates['X'] = df_Elastic_FinishDay['X']
        df_Elastic_DiffBetweenTwoDates['Y'] = df_Elastic_FinishDay['Y']
        df_Elastic_DiffBetweenTwoDates['Z'] = df_Elastic_FinishDay['Z']
        # Add columns for relative displacements in between the two dates
        df_Elastic_DiffBetweenTwoDates['RelDispX'] = df_Elastic_FinishDay['DispX'] - df_Elastic_StartDay['DispX']
        df_Elastic_DiffBetweenTwoDates['RelDispY'] = df_Elastic_FinishDay['DispY'] - df_Elastic_StartDay['DispY']
        df_Elastic_DiffBetweenTwoDates['RelDispZ'] = df_Elastic_FinishDay['DispZ'] - df_Elastic_StartDay['DispZ']
    else:
        raise ValueError("The first 3 columns are not identical between the two Elastic DataFrames.")

    ##############################################################################
    ###load and process data for viscous simulation (one transient simulation)####
    ##############################################################################
    ####Find corresponding day number since beginning of simu for viscous case (one transient simulations)
    day_start = int(StartDate_str[:2])
    month_start = int(StartDate_str[2:4])
    year_start = int(StartDate_str[4:])
    StartDate = date(year_start, month_start, day_start)
    StartDayNumber = (StartDate - RefStartDate).days ##That's the start day number knowing that 01/07/2011 is day 1

    day_finish = int(FinishDate_str[:2])
    month_finish = int(FinishDate_str[2:4])
    year_finish = int(FinishDate_str[4:])
    FinishDate = date(year_finish, month_finish, day_finish)
    FinishDayNumber = (FinishDate - RefStartDate).days ##That's the finish day number knowing that 01/07/2011 is day 1

    ####Open file
    if FinishDayNumber < -255: ###This is simu called Pumping 2010
        StartSimuDayNumber_ViscousSimu = -315  ###That's the day number corresponding to the beginning of simu pumping 2010
        filename_viscous = './ScalarOutput/SurfaceOutput_Prono_NoD_PCav_P2010_tsp1d_Out5d_Tmap_.dat'
        df_Viscous = pd.read_csv(Pathroot.joinpath(filename_viscous), names=Col_Names_Viscous, delim_whitespace=True)
        df_Viscous['TimestepSize'] = 1 ##This simu as a timestep of one day
        df_Viscous['DayofSimu'] = df_Viscous['TimestepSize']* (df_Viscous['Timestep']-1 )+StartSimuDayNumber_ViscousSimu
    elif -255 < StartDayNumber < 85 and FinishDayNumber < 85: ###This is simu called Refill20102011
        StartSimuDayNumber_ViscousSimu = -255  ###That's the day number corresponding to the beginning of simu refill 20102011
        filename_viscous = './ScalarOutput/SurfaceOutput_Prono_NoD_PCav_Refill20102011_tsp5d_Out30d_Tmap_.dat'
        df_Viscous = pd.read_csv(Pathroot.joinpath(filename_viscous), names=Col_Names_Viscous, delim_whitespace=True)
        df_Viscous['TimestepSize'] = 5 ##This simu as a timestep of 5 days
        df_Viscous['DayofSimu'] = df_Viscous['TimestepSize'] * (df_Viscous['Timestep'] - 1) + StartSimuDayNumber_ViscousSimu
    elif StartDayNumber < 85 and 85 < FinishDayNumber < 105: ##Startday is in simu Refill20102011, Finish day is in simu Pumping2011
        StartSimuDayNumber_ViscousSimu_start = -255 ###That's the day number corresponding to the beginning of simu refill 20102011
        filename_viscous_start = './ScalarOutput/SurfaceOutput_Prono_NoD_PCav_Refill20102011_tsp5d_Out30d_Tmap_.dat'
        df_Viscous_start = pd.read_csv(Pathroot.joinpath(filename_viscous_start), names=Col_Names_Viscous, delim_whitespace=True)
        df_Viscous_start['TimestepSize'] = 5 ##This simu as a timestep of 5 days
        df_Viscous_start['DayofSimu'] = df_Viscous_start['TimestepSize'] * (df_Viscous_start['Timestep'] - 1) + StartSimuDayNumber_ViscousSimu_start
        StartSimuDayNumber_ViscousSimu_finish = 85 ###That's the day number corresponding to the beginning of simu pumping 2011
        filename_viscous_finish = './ScalarOutput/SurfaceOutput_Prono_NoD_PCav_P2011_tsp1d_Out5d_Tmap_.dat'
        df_Viscous_finish = pd.read_csv(Pathroot.joinpath(filename_viscous_finish), names=Col_Names_Viscous, delim_whitespace=True)
        df_Viscous_finish['TimestepSize'] = 1 ##This simu as a timestep of one day
        df_Viscous_finish['DayofSimu'] = df_Viscous_finish['TimestepSize'] * (df_Viscous_finish['Timestep'] - 1) + StartSimuDayNumber_ViscousSimu_finish
        df_Viscous = pd.concat([df_Viscous_start, df_Viscous_finish], ignore_index=True)
    elif 85 <StartDayNumber < 105 and 85 < FinishDayNumber < 105:###This is simu called Pumping2011
        StartSimuDayNumber_ViscousSimu = 85###That's the day number corresponding to the beginning of simu pumping 2011
        filename_viscous = './ScalarOutput/SurfaceOutput_Prono_NoD_PCav_P2011_tsp1d_Out5d_Tmap_.dat'
        df_Viscous = pd.read_csv(Pathroot.joinpath(filename_viscous), names=Col_Names_Viscous, delim_whitespace=True)
        df_Viscous['TimestepSize'] = 1##This simu as a timestep of one day
        df_Viscous['DayofSimu'] = df_Viscous['TimestepSize']* (df_Viscous['Timestep']-1 )+StartSimuDayNumber_ViscousSimu
    elif 85 < StartDayNumber < 105  and 105 < FinishDayNumber < 450: ##Startday is in simu Pumping2011, Finish day is in simu Refill20112012
        StartSimuDayNumber_ViscousSimu_start = 85 ###That's the day number corresponding to the beginning of simu Pumping 011
        filename_viscous_start = './ScalarOutput/SurfaceOutput_Prono_NoD_PCav_P2011_tsp1d_Out5d_Tmap_.dat'
        df_Viscous_start = pd.read_csv(Pathroot.joinpath(filename_viscous_start), names=Col_Names_Viscous, delim_whitespace=True)
        df_Viscous_start['TimestepSize'] = 1 ##This simu as a timestep of one day
        df_Viscous_start['DayofSimu'] = df_Viscous_start['TimestepSize'] * (df_Viscous_start['Timestep'] - 1) + StartSimuDayNumber_ViscousSimu_start
        StartSimuDayNumber_ViscousSimu_finish = 105 ###That's the day number corresponding to the beginning of simu refill20112012
        filename_viscous_finish = './ScalarOutput/SurfaceOutput_Prono_NoD_PCav_Refill20112012_tsp5d_Out30d_Tmap_.dat'
        df_Viscous_finish = pd.read_csv(Pathroot.joinpath(filename_viscous_finish), names=Col_Names_Viscous, delim_whitespace=True)
        df_Viscous_finish['TimestepSize'] = 5 ##This simu as a timestep of 5 days
        df_Viscous_finish['DayofSimu'] = df_Viscous_finish['TimestepSize'] * (df_Viscous_finish['Timestep'] - 1) + StartSimuDayNumber_ViscousSimu_finish
        df_Viscous = pd.concat([df_Viscous_start, df_Viscous_finish], ignore_index=True)
    else:
        raise ValueError("The chosen starting and/or finishing dates are not treated in this code. More developments required")
    ###To be safe we re-initialize the indexes of df_Viscous
    df_Viscous.reset_index(drop=True)
    ###We need to find the starting and finishing dates in the simu outputs (if days are fitting)
    SimuDays = df_Viscous['DayofSimu'].drop_duplicates()
    ###Check whether starting and finishing dates are among simu days and if not interpolate velocities at these days
    for obsday in [StartDayNumber, FinishDayNumber]:
        if obsday in SimuDays.values:##If start day or finish day are already in the simu days, this operation is not required
            continue
        print('day= ',obsday,'is not in simulation days. Proceeding to temporal interpolation.')
        ### find days in simu surrounding the day of interest
        simuday_before = df_Viscous[df_Viscous['DayofSimu'] <obsday]['DayofSimu'].iloc[-1]
        simuday_after = df_Viscous[df_Viscous['DayofSimu'] >obsday]['DayofSimu'].iloc[0]
        ###Check that nodes are ordered the same for both days
        df_Viscous_before = df_Viscous[df_Viscous['DayofSimu'] == simuday_before].reset_index(drop=True)
        df_Viscous_after = df_Viscous[df_Viscous['DayofSimu'] == simuday_after].reset_index(drop=True)
        if df_Viscous_before['NodeNumber'].equals(df_Viscous_after['NodeNumber']):
            print("The 'NodeNumber' values are identical for the two days.")
        else:
            raise ValueError("The 'NodeNumber' values differ between the two days.")
        # Copy the first four columns from df_Viscous_before [tsp, tsp2,
        df_Viscous_weighted = df_Viscous_before.iloc[:, :4].copy()
        # Compute the weighted average between day before and day after for the remaining columns
        for col in df_Viscous_before.columns[4:]:
            weight = (obsday-simuday_before)/(simuday_after-simuday_before)
            df_Viscous_weighted[col] = ((1-weight) * df_Viscous_before[col] + weight * df_Viscous_after[col])
        ### fill-up df_Viscous with this new df for the day of interrest at the right position
        last_row_index = df_Viscous[df_Viscous['DayofSimu'] == simuday_before].index[-1]
        # insert the new df
        top_df_Viscous = df_Viscous.iloc[:last_row_index + 1]  # Partie avant l'insertion
        bottom_df_Viscous = df_Viscous.iloc[last_row_index + 1:]  # Partie après l'insertion
        df_Viscous = pd.concat([top_df_Viscous, df_Viscous_weighted, bottom_df_Viscous], ignore_index=True)
        df_Viscous.reset_index(drop=True)
    ###To get the total displacement in between the two dates one needs to integrate the velocity field !!!
    # The idea is to initialize particles at mesh node and then to track displacements of particles on the considered period
    # Initialize Particles Coords:one line per tsp, one column per particle (i.e. node), for x, y and z
    NewSimuDays = df_Viscous['DayofSimu'].drop_duplicates().reset_index(drop=True)
    NumberOfSimuDays = NewSimuDays[NewSimuDays == FinishDayNumber].index[0]-NewSimuDays[NewSimuDays == StartDayNumber].index[0]
    Particles_Coords = np.zeros((NumberOfSimuDays+1, len(df_Viscous[df_Viscous['DayofSimu'] == StartDayNumber]), 3)) # 3 components (x, y, z) for each point
    # For first timestep of considered period, particles coordinates are coordinate of mesh nodes:
    filtered_df = df_Viscous[df_Viscous['DayofSimu'] == StartDayNumber]
    Particles_Coords_Start = filtered_df[['X', 'Y', 'Z']].values
    Particles_Coords[0, :, :] = Particles_Coords_Start
    ###TO CHECK FROM HERE
    ###Now loop on time steps of the considered period (one tsp per day) to displace particle day after day
    for k,day_idx in enumerate(range(NewSimuDays[NewSimuDays == StartDayNumber].index[0],NewSimuDays[NewSimuDays == FinishDayNumber].index[0])):
        k =k+1##as k=0 is for initial condition
        day_prev = NewSimuDays.loc[day_idx]
        day_curr = NewSimuDays.loc[day_idx+1]
        # Get simu output data for timesteps of interest
        v_prev = df_Viscous[df_Viscous['DayofSimu']==day_prev]
        v_curr = df_Viscous[df_Viscous['DayofSimu']==day_curr]
        # Extract velocity components
        velocity_prev = v_prev[['U', 'V', 'W']].values
        velocity_curr = v_curr[['U', 'V', 'W']].values
        # Get previous positions of particles
        Particles_Coords_prev = Particles_Coords[k-1, :, :]
        # Interpolate velocity fields at particles previous location
        v_prev_interp = griddata(Particles_Coords_Start, velocity_prev, Particles_Coords_prev, method='linear')
        v_curr_interp = griddata(Particles_Coords_Start, velocity_curr, Particles_Coords_prev, method='linear')
        # Compute new positions of particles using the trapezoidal rule
        Particles_Coords[k, :, :] = Particles_Coords[k-1, :, :] + 0.5 * (v_prev_interp + v_curr_interp)*((day_curr-day_prev)/365.25)
    ### Compute total displacement field as the difference between particle positions at last tsp and initial particle positions
    Total_Displacement_Viscous=Particles_Coords[-1, :, :]-Particles_Coords_Start
    ### Build panda dataframe to store displacements together with all required information
    df_Viscous_DiffBetweenTwoDates = df_Viscous[df_Viscous['DayofSimu'] == StartDayNumber].iloc[:, :7].copy()
    # Add columns for relative displacements in between the two dates
    df_Viscous_DiffBetweenTwoDates['RelDispX'] = Total_Displacement_Viscous[:, 0]
    df_Viscous_DiffBetweenTwoDates['RelDispY'] = Total_Displacement_Viscous[:, 1]
    df_Viscous_DiffBetweenTwoDates['RelDispZ'] = Total_Displacement_Viscous[:, 2]
    ###Remove lines with NaN due to interpolation (it's always nodes that are not above cavity so not a big problem)
    df_Viscous_DiffBetweenTwoDates=df_Viscous_DiffBetweenTwoDates.dropna()

    ##############################################
    # NOW START THE FIGURE STRICTLY SPEAKING #####
    ##############################################
    for j,case in enumerate(['Elastic', 'Viscous', 'Obs']): ##loop on case (one per column)
        ### get corresponding data frame
        if case == 'Elastic':
            df = df_Elastic_DiffBetweenTwoDates
        elif case == 'Viscous':
            df = df_Viscous_DiffBetweenTwoDates
        elif case == 'Obs':
            df = df_Obs_DiffBetweenTwoDates
        ### loop on the three dimensions x,y,z
        for i,dim in enumerate(['X','Y','Z']):
            print('Case =', case, 'for', dim)
            ###get proper dimension
            displacementfield_name='RelDisp{}'.format(dim)
            ###interpolate displacement of interest on the refined grid
            Interpolated_disp = Interpolate_field(df,displacementfield_name,X,Y)
            ###get proper ax
            ax = axes[i,j]
            ax.set_xlim([0, xmax-xmin])
            ax.set_ylim([0, ymax-ymin])
            if i==0 and j==0:
                ax.set_title('Model Elastic', fontsize=21, weight='bold')
            elif i==0 and j==1:
                ax.set_title('Model Viscous', fontsize=21, weight='bold')
            elif i==0 and j==2:
                ax.set_title('Stake measurements', fontsize=21, weight='bold')

            #######################################
            ##########    START PLOT    ###########
            #######################################
            if i ==0 : ##Ux
                if j == 0: ##elastic
                    clevs = np.arange(dx_min_el, dx_max_el, lev_ticks_range_dx_el/color_resolution)  ## cbar for shading
                    levs_ticks = np.arange(dx_min_el, dx_max_el, lev_ticks_range_dx_el)
                    norm = norm_dx_el
                else: ##same scale for visc and obs
                    clevs = np.arange(dx_min_obs, dx_max_obs, lev_ticks_range_dx_obs/color_resolution)  ## cbar for shading
                    levs_ticks = np.arange(dx_min_obs, dx_max_obs, lev_ticks_range_dx_obs)
                    norm = norm_dx_obs
            elif i == 1: ##Uy
                if j == 0: ##elastic
                    clevs = np.arange(dy_min_el, dy_max_el, lev_ticks_range_dy_el/color_resolution)  ## cbar for shading
                    levs_ticks = np.arange(dy_min_el, dy_max_el, lev_ticks_range_dy_el)
                    norm = norm_dy_el
                else: ##same scale for visc and obs
                    clevs = np.arange(dy_min_obs, dy_max_obs, lev_ticks_range_dy_obs/color_resolution)  ## cbar for shading
                    levs_ticks = np.arange(dy_min_obs, dy_max_obs, lev_ticks_range_dy_obs)
                    norm = norm_dy_obs
            elif i ==2: ##Uz
                if j == 0: ##elastic
                    clevs = np.arange(dz_min_el, dz_max_el, lev_ticks_range_dz_el/color_resolution)  ## cbar for shading
                    levs_ticks = np.arange(dz_min_el, dz_max_el, lev_ticks_range_dz_el)
                    norm = norm_dz_el
                else: ##same scale for visc and obs
                    clevs = np.arange(dz_min_obs, dz_max_obs, lev_ticks_range_dz_obs/color_resolution)  ## cbar for shading
                    levs_ticks = np.arange(dz_min_obs, dz_max_obs, lev_ticks_range_dz_obs)
                    norm = norm_dz_obs


            ###Fills up the map with colors for displacement field interpolated
            CS1 = ax.contourf(X-xmin, Y-ymin, Interpolated_disp*100, clevs, cmap=cmap, extend='both', norm=norm)
            if j == 2: ###For observations only:
                ###Plot stake positions
                ax.scatter(df_Obs_DiffBetweenTwoDates['X'] -xmin, df_Obs_DiffBetweenTwoDates['Y'] -ymin, color='k', marker='.', linewidths=0.4)
                # Add labels to each scatter point
                for l, label in enumerate(df_Obs_DiffBetweenTwoDates['Stake']):  # Iterate through each point and label
                    ax.text(df_Obs_DiffBetweenTwoDates['X'].iloc[l] -xmin, df_Obs_DiffBetweenTwoDates['Y'].iloc[l] -ymin,label,fontsize=9,ha='right',va='bottom')
            ###Show colorbar only for first and last column
            if j == 0 or j == 2: ###For observations only:
                cbar = plt.colorbar(CS1, ax=ax, ticks=levs_ticks, orientation='vertical', label=r'$U_{}$ [cm]'.format(dim))

            ###Plot cavity contour
            ax.plot(cavity_contour[:, 0]-xmin,cavity_contour[:, 1]-ymin,color='k',linestyle='-',linewidth=0.9)
            ###Plot also the one from sonar
            ax.plot(Df_GLsonar['X'].values-xmin,Df_GLsonar['Y'].values-ymin,color='gray',linestyle='-',linewidth=1.3)
    ################################################################################
    # SAVE THE FIGURES #####
    ################################################################################
    # name_output_fig = 'RelDensProfile_Forage{}_{}_DInitHL_'.format(ForageNumber,  Case)
    # name_path = '/home/brondexj/BETTIK/MyTaconnaz/MyTaco_Porous_DensSemiLagVSModifHeatSolv/Figures/Case_FullSimu_TestNbIntTsp/'
    # path_output_fig = Path(name_path)
    # fig.savefig(path_output_fig.joinpath(name_output_fig))

    plt.show()