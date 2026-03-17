"""
@author: jbrondex

Description:
------------
This file aims at plotting:
 1/ four subplots corresponding to the maps of damage for combination 7
 at the end of 2010 pumping for the failure criterion Max Princip. Stress, Hayurst, von Mises Coulomb

"""

################################################################################
# import #####
################################################################################
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cmx
from matplotlib.path import Path as mpltPath
from matplotlib.patches import PathPatch

from pathlib import Path
from datetime import datetime, date, timedelta ###To handle dates
from matplotlib.dates import MonthLocator, DayLocator, HourLocator, DateFormatter, drange
import numpy as np
from scipy.interpolate import griddata
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
#####################################################################################
# Define function to interpolate field defined on mesh nodes to the line profile ####
#####################################################################################
def Interpolate_field(df, field_name, x, y): ###Returned field interpolated from mesh grid to point (x,y) of transect
    # x, y, and z coordinates of the dataframe
    xi = df['X'].values
    yi = df['Y'].values
    field = df[field_name].values
    # Interpolate the value of field
    result = griddata((xi, yi), field, (x, y), rescale=True)
    return result

def Coord_transect(pt1, pt2): ###return list of coordinates corresponding to line passing by pt1 and pt2
    x1, y1 = pt1 ##pt1 defining transect
    x2, y2 = pt2 ##pt2 defining transect
    xmin=min(x1, x2)
    xmax=max(x1, x2)
    x = np.linspace(xmin, xmax, 1000)
    y = ((y2-y1)/(x2-x1))*(x-x1)+y1
    coord = [x,y]
    return coord

def Distance_along_transect(pt1, pt2, pt3): ###return distance of pt3 along transect (pt1,pt2)
    x1, y1 = pt1 ##pt1 defining transect
    x2, y2 =pt2 ##pt2 defining transect
    x, y = pt3
    x0=min(x1, x2)
    if x0 == x1:
        y0 = y1
    elif x0 == x2:
        y0 = y2
    dist_along_line = np.sqrt((x-x0)**2+(y-y0)**2)
    return dist_along_line

def Distance_to_line(pt1, pt2, pt3): ###return distance of pt3 to line define by points (pt1, pt2)
    x1, y1 = pt1 ##pt1 defining transect
    x2, y2 = pt2 ##pt2 defining transect
    x3, y3 = pt3 ##pt3 out of transect. How far is it from transect ?

    det = (x2-x1)*(y3-y1)-(y2-y1)*(x3-x1) ##This is det(AB;AC) with A and B the two extremities of transect and C is pt3
    distance = det/np.sqrt((x2-x1)**2 + (y2-y1)**2) ##det(AB;AC)/||AB|| = ||AC||*sin(AB,AC) = ||CD||
    return distance

def Proj_on_line(pt1, pt2, pt3):##Return coords of pt3 projected on line defined by (pt1,pt2)
    x1, y1 = pt1
    x2, y2 = pt2
    x3, y3 = pt3
    dx, dy = x2-x1, y2-y1
    scalarprod = dy*(y3-y1)+dx*(x3-x1) ##AB.AC=||AB||*||AC||*cos(AB,AC)=||AB||*||AD||
    a = scalarprod/(dx*dx + dy*dy) ## a = ||AB||*||AD|| / ||AB||**2  = ||AD|| / ||AB||
    return x1+a*dx, y1+a*dy

##Function below is used to sort point of contour clockwise
def sort_clockwise(list_of_xy_coords):
    cx, cy = list_of_xy_coords.mean(0)
    x, y = list_of_xy_coords.T
    angles = np.arctan2(x-cx, y-cy)
    indices = np.argsort(angles)
    return list_of_xy_coords[indices]

################################################################################
# prepare plots #####
################################################################################
# global figure pimping
# plt.rc('text', usetex=True)
plt.rc('xtick', labelsize=22)
plt.rc('ytick', labelsize=22)
plt.rc('axes', labelsize=22)
plt.rc('legend', fontsize=26)


if __name__ == "__main__":
    ####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####
    ####~~~~~~~~~~~   LOAD ALL DATA AND PLOT DIRECTLY IN THIS PART OF THE CODE   ~~~~~~~~~~~~~####
    ####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####
    ### General parameter of simu
    RefStartDate = date(int(2011), 7, 1) ##Every days are numbered relative to 1st July 2011 (day 1)
    StartDayNumber_Pumping2010 = "-315"
    StartYear_Pumping2010 = "2010"
    StartDate_Pumping2010 = RefStartDate + timedelta(days=int(StartDayNumber_Pumping2010))
    ###Provide coordinate of points defining lines corresponding to transect on which we want to plot SigmaI: one couple of point per transect
    List_Transect=["AA'","BB'", "CC'"]
    List_Coord_pt1=[[947954.6,2105001.5], [947976.1, 2105120.6], [947943.3, 2105052.1]]
    List_Coord_pt2=[[948027.9,2105114.9], [948024.0,2104971.5], [948112.8,2105049.4]]
    ### Step in the full process from 2010 to 2014
    Step_List = ['P2010', 'Refill20102011']  # , 'P2011', 'Refill20112012', 'P2012', 'Refill20122013']
    StepTsp_List = [1, 5]  # , 1, 5, 1, 5] ##Timestep size (in days) corresponding to simulation step (Step_List)
    StepOut_List = [5, 30]  # , 5, 30, 5, 30] ##Output interval (in days) corresponding to simulation step (Step_List)
    ## Where pressure applies: cavity only: 'PCavityOnly', restricted to a conduit: 'PRestricted', everywhere above cold/temperate transition: 'PNotRestric'
    Case = 'PCav'
    T = 'Tmap' ##'Temperate', 'Tminus2' or 'Tmap'
    Col_Crevasses_Circ = '#00ff00' ###'#FF0000' ###'#800080' ###Color for representation of circular crevasses
    Col_Crevasses_Other = '#800080' ###Color for representation of non circular crevasses
    Col_Cavity = '#00FFFF'###color for representation of cavity
    Col_Transect = 'w'###color for representation of transect

    Colormap_for_stressmap = 'magma' ###Colorm map to choose for map of stress
    cmap = cmx.get_cmap(Colormap_for_stressmap )
    #cmap = colors.LinearSegmentedColormap.from_list("viridis_darker", cmap(np.linspace(0.2, 1, 256)))

    Colormap_for_Simudays = 'copper_r'
    # cmap_simudays = cmx.get_cmap(Colormap_for_Simudays )
    # cmap_simudays = colors.LinearSegmentedColormap.from_list("cividis_darker", cmap_simudays(np.linspace(0.2, 1, 256)))
    ################################################################
    ####  OPEN DATA CORRESPONDING TO CREVASSES AND GROUNDEDMASK ####
    ################################################################
    Pathroot_Obs = Path('/home/brondexj/BETTIK/TeteRousse/MyTeterousse_GeoGag/Data/.')
    Filename_Crevasses = 'crevasses_xyz_numbered.dat'
    ###load file as dataframe
    Df_Crevasses = pd.read_csv(Pathroot_Obs.joinpath(Filename_Crevasses), names=['X', 'Y', 'Z', 'Crevasse Number', 'IsCircular'], delim_whitespace=True)
    ###Also load glacier contour
    Filename_GlacierContour = 'Contour_TeteRousse2012.dat'
    Df_GlacierContour = pd.read_csv(Pathroot_Obs.joinpath(Filename_GlacierContour), names=['X', 'Y'], delim_whitespace=True)
    xc = Df_GlacierContour['X'].values
    xc=np.append(xc,xc[0]) ##step required to close the contour
    yc = Df_GlacierContour['Y'].values
    yc = np.append(yc, yc[0])  ##step required to close the contour
    ###Open Data corresponding to grounded mask
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
    #########################################
    ####  OPEN OUTPUT OF THE SIMULATIONS:####
    #########################################
    ### Load files containing surface output of the simulations:
    Pathroot_SimuOutput = Path('/home/brondexj/BETTIK/TeteRousse/MyTeterousse_GeoGag/ScalarOutput/.')
    ###NAMES OF OUTPUT DATA (same order as in dat.names file)
    ###BE WARE: I use the convention SigmaI>SigmaII>SigmaIII which is not the convention of the EigenStress solver in Elmer
    Col_Names_NoD = ['Timestep', 'DayOfSimu', 'BC', 'NodeNumber', 'X', 'Y', 'Z', 'U', 'V', 'W','SigmaI', 'D','Chi']
    ###~~~~~~~Prepare the figure 1 : subplot with the surface eq stress for the four criterion tested
    ###Prepare the subplot
    nrows = 2  ###Sigma_I, Uz
    ncols = 2  ###1GPa, 9GPa, Diff
    fig1, axes = plt.subplots(nrows, ncols, figsize=(15, 12), sharex=True, sharey=True) #, constrained_layout=True)  # , gridspec_kw={'hspace': -0.15})
    plt.subplots_adjust(bottom=0.2,wspace=0.1)
    ### Start loop on the various criterion considered
    for idx,(Criterion,Title) in enumerate(zip(['MPS', 'vM', 'Ha', 'Co'], ['Max. principal stress', 'von Mises', 'Hayhurst', 'Coulomb'])):
        ###We create a single dataframe for all considered steps of the simu
        Data_Simu_NoD = pd.DataFrame()  ##Initialize an empty dataframe
        for i, (Step, StepTsp, StepOut) in enumerate(zip(Step_List, StepTsp_List, StepOut_List)):
            filename = 'SurfaceOutput_Prono_B172_Sth01_Lh00_{}_{}_tsp{}d_Out{}d_Crit{}_{}_.dat'.format(Case, Step, str(StepTsp), str(StepOut), Criterion, T)
            print('Opening file:', filename)
            Data_Simu_NoD_tmp = pd.read_csv(Pathroot_SimuOutput.joinpath(filename), names=Col_Names_NoD, delim_whitespace=True)
            ###Drop duplicate lines (bug of SaveLine solver)
            Data_Simu_NoD_tmp.drop_duplicates(inplace=True)
            ###Set Day of Simu counting from begining of simu corresponding to step Pump2010 (Step 1)
            if i == 0:
                Data_Simu_NoD_tmp['DayOfSimu'] = Data_Simu_NoD_tmp['DayOfSimu'] * StepTsp
            else:
                Data_Simu_NoD_tmp['DayOfSimu'] = np.max(Data_Simu_NoD['DayOfSimu']) + Data_Simu_NoD_tmp['DayOfSimu'] * StepTsp
            ###We create an additionnal column storing the step of simu we are dealing with
            Data_Simu_NoD_tmp['Step'] = Step
            data = [Data_Simu_NoD, Data_Simu_NoD_tmp]
            Data_Simu_NoD = pd.concat(data, ignore_index=True)
        ### Get the last day of the step pumping 2010
        Df_LastDayPump2010=Data_Simu_NoD[Data_Simu_NoD['DayOfSimu']==np.max(Data_Simu_NoD[Data_Simu_NoD['Step']=='Refill20102011']['DayOfSimu'])].copy()
        ########################################################################################
        ###     PLOT MAP OF DAMAGE AT THE SURFACE AT THE END OF PUMPING FOR EACH CRITERION   ###
        ########################################################################################
        ### Get proper ax
        row, col = divmod(idx, ncols)  # Get row, col index from the idx
        ax = axes[row, col]  # Select the axis in the grid
        if row == 1:
            ax.set_xlabel(r'X [km]', fontsize=22)
        if col == 0:
            ax.set_ylabel(r'Y [km]', fontsize=22)
        ax.tick_params(labelsize=18)  # fontsize of the tick labels
        # ax.set_xlabel("")
        # ax.set_ylabel("")
        # ax.set_xticklabels([])  # Remove x labels
        # ax.set_yticklabels([])  # Remove y labels
        ax.grid(True)
        ax.grid(alpha=0.5)
        ###Force x and y axis to same scale
        ax.set_aspect('equal', adjustable='box')
        ax.set_title(Title, fontsize=21, weight='bold')
        # plt.xlim([(np.floor(np.min(Df_LastDayPump2010['X']))-10)/1000, (np.floor(np.max(Df_LastDayPump2010['X']))+10)/1000])
        # plt.ylim([(np.floor(np.min(Df_LastDayPump2010['Y']))-10)/1000, (np.floor(np.max(Df_LastDayPump2010['Y']))+10)/1000])
        ###Create refined regular grid and interpolate field over grid
        x = np.arange(np.floor(np.min(Df_LastDayPump2010['X']))-10,np.floor(np.max(Df_LastDayPump2010['X']))+10,0.1)
        y = np.arange(np.floor(np.min(Df_LastDayPump2010['Y']))-10,np.floor(np.max(Df_LastDayPump2010['Y']))+10,0.1)
        X, Y = np.meshgrid(x, y)
        Damage = Interpolate_field(Df_LastDayPump2010,'D',X,Y)
        #shading
        clevs=np.arange(0.0,0.10001,0.01) ## cbar for shading
        cmap=Colormap_for_stressmap
        #colorbar
        levs_ticks=np.arange(0.0,0.1001,0.025)
        ###Fills up the map with colors for SigmaEq
        CS1 = ax.contourf(X/1000, Y/1000, Damage, clevs, cmap=cmap,extend='both')
        ax.plot(xc / 1000, yc / 1000, color='k', linewidth=2)
        ###Show colorbar
        # cbar = ax.colorbar(CS1, ticks=levs_ticks, orientation='vertical', label=r'$\sigma_\mathrm{eq}$ [kPa]')
        ###Below we remove colors that are outside of glacier contour
        clippath = mpltPath(np.c_[xc/1000, yc/1000])
        patch = PathPatch(clippath, facecolor='none')
        ax.add_patch(patch)
        for c in [CS1]:
            c.set_clip_path(patch)
        ###Plot cavity contour
        ax.plot(cavity_contour[:, 0]/1000,cavity_contour[:, 1]/1000,color=Col_Cavity,linestyle='-',linewidth=2.6)
        ###plot crevasses as continuous line
        for crev_num in np.arange(Df_Crevasses['Crevasse Number'].min(), Df_Crevasses['Crevasse Number'].max() + 1):
            ###get proper points
            Df_plot = Df_Crevasses[Df_Crevasses['Crevasse Number'] == crev_num]
            if Df_plot['IsCircular'].all():  ##different colors for circular crevasses and other crevasses
                col = Col_Crevasses_Circ
            else:
                continue #col = Col_Crevasses_Other
            ax.plot(Df_plot['X'].values/1000, Df_plot['Y'].values/1000, color=col, linestyle='-', linewidth=2)
    # Add a common colorbar to all subplots below the figure
    cbar_ax = fig1.add_axes([0.15, 0.082, 0.7, 0.027])  # [left, bottom, width, height]
    fig1.colorbar(CS1, ticks=levs_ticks, cax=cbar_ax, orientation='horizontal', label=r'Damage')
    ###Show map
    plt.show()





