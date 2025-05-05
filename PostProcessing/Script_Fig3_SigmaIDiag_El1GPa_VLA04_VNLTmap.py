"""
@author: jbrondex

Description:
------------
This file aims at plotting the surface maximum principal stress obtained with the elastic (E=1GPa),
viscous line (A=0.4 MPa-1 a-1), and viscous non line (T = Tmap) frameworks. Each time the background sigma_I
(i.e. stress obtained when there is no cavity) is substracted to the stress field with cavity."""

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
    Col_Crevasses_Circ = '#FF0000' ###'#800080' ###Color for representation of circular crevasses
    Col_Crevasses_Other = '#800080' ###Color for representation of non circular crevasses
    Col_Cavity = '#00FFFF'###color for representation of cavity
    Col_Transect = 'lightgrey'###color for representation of transect

    Colormap_for_stressmap = 'viridis' ###Colorm map to choose for map of stress
    cmap = cmx.get_cmap(Colormap_for_stressmap )
    cmap = colors.LinearSegmentedColormap.from_list("viridis_darker", cmap(np.linspace(0.2, 1, 256)))
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
    ###BE WARE: I use the convention SigmaI>SigmaII>SigmaIII which is not the convention of the EigenStress solver in Elmer viscous stress
    Col_Names_Viscous = ['Timestep', 'BC', 'NodeNumber', 'X', 'Y', 'Z', 'U', 'V', 'W', 'Pressure', 'SigmaIII','SigmaII','SigmaI', 'temperature', 'Sxx', 'Syy', 'Szz', 'Sxy', 'Sxz', 'Syz']
    Col_Names_Elastic = ['Timestep', 'BC', 'NodeNumber', 'X', 'Y', 'Z', 'Dispx', 'DispY', 'DispZ', 'vonMises', 'SigmaI','SigmaII','SigmaIII']
    ###We create a single dataframe for each constitutive law : one for elastic, one for viscous line, one for viscous non line
    Data_Simu_El = pd.DataFrame()  ##Initialize an empty dataframe
    Data_Simu_VL = pd.DataFrame()
    Data_Simu_VNL = pd.DataFrame()
    for case in ['NoCavity', 'EmptyCavity']:
        filename_el = 'SurfaceOutput_DiagElastic_Poisson03_E1GPa_{}_NoDispLat_.dat'.format(case)
        filename_vl = 'SurfaceOutput_DiagStokes_Line_{}_A04_NoSliding_NoDispLat.dat'.format(case)
        filename_vnl = 'SurfaceOutput_DiagStokes_NonLine_{}_Tmap_NoSliding_NoDispLat.dat'.format(case)

        print('Opening file:', filename_el)
        Data_Simu_El_tmp = pd.read_csv(Pathroot_SimuOutput.joinpath(filename_el), names=Col_Names_Elastic, delim_whitespace=True)
        Data_Simu_El_tmp['Case'] = case
        Data_Simu_El = pd.concat([Data_Simu_El, Data_Simu_El_tmp], ignore_index=True)

        print('Opening file:', filename_vl)
        Data_Simu_VL_tmp = pd.read_csv(Pathroot_SimuOutput.joinpath(filename_vl), names=Col_Names_Viscous, delim_whitespace=True)
        Data_Simu_VL_tmp['Case'] = case
        Data_Simu_VL = pd.concat([Data_Simu_VL, Data_Simu_VL_tmp], ignore_index=True)

        print('Opening file:', filename_vnl)
        Data_Simu_VNL_tmp = pd.read_csv(Pathroot_SimuOutput.joinpath(filename_vnl), names=Col_Names_Viscous, delim_whitespace=True)
        Data_Simu_VNL_tmp['Case'] = case
        Data_Simu_VNL = pd.concat([Data_Simu_VNL, Data_Simu_VNL_tmp], ignore_index=True)

    # #####################################################################
    # ###     PLOT MAP OF FIRST MAXI PRINCIPAL STRESS AT THE SURFACE    ###
    # #####################################################################
    ###Create refined regular grid on which all SigmaI for all cases will be interpolated
    x = np.arange(np.floor(np.min(Data_Simu_El['X'])) - 10,np.floor(np.max(Data_Simu_El['X'])) + 10, 0.1)
    y = np.arange(np.floor(np.min(Data_Simu_El['Y'])) - 10, np.floor(np.max(Data_Simu_El['Y'])) + 10, 0.1)
    X, Y = np.meshgrid(x, y)
    ###Prepare the subplot
    nrows = 1
    ncols = 3  ###elastic, viscous, observed
    fig, axes = plt.subplots(nrows, ncols, figsize=(15, 12), sharex=True, sharey=True, constrained_layout=True)

    ###### START PLOT ######
    for k,law in enumerate(['El','VL','VNL']):
        df_cav =pd.DataFrame()
        df_nocav=pd.DataFrame()
        if law == 'El':
            df_nocav =  Data_Simu_El[Data_Simu_El['Case']=='NoCavity']
            df_cav = Data_Simu_El[Data_Simu_El['Case']=='EmptyCavity']
            Title = 'Elastic'
        elif law == 'VL':
            df_nocav =  Data_Simu_VL[Data_Simu_VL['Case']=='NoCavity']
            df_cav = Data_Simu_VL[Data_Simu_VL['Case']=='EmptyCavity']
            Title = 'Linear Viscous'
        elif law=='VNL':
            df_nocav =  Data_Simu_VNL[Data_Simu_VNL['Case']=='NoCavity']
            df_cav = Data_Simu_VNL[Data_Simu_VNL['Case']=='EmptyCavity']
            Title = 'Glen\'s law'
        ##get proper ax
        ax = axes[k]
        ax.set_xlabel(r'X [km]', fontsize=22)
        ax.set_ylabel(r'Y [km]', fontsize=22)
        ax.tick_params(labelsize=18)  # fontsize of the tick labels
        ax.set_xlabel("")
        ax.set_ylabel("")
        ax.set_xticklabels([])  # Remove x labels
        ax.set_yticklabels([])  # Remove y labels
        ax.grid(True)
        ax.grid(alpha=0.5)
        ax.set_title(Title, fontsize=21, weight='bold')
        ###Force x and y axis to same scale
        ax.set_aspect('equal', adjustable='box')
        ###Poject SigmaI on interpolation grid
        SigmaI_NoCav = Interpolate_field(df_nocav, 'SigmaI', X, Y)
        SigmaI_Cav = Interpolate_field(df_cav, 'SigmaI', X, Y)
        ###shading
        clevs=np.arange(-0.01,121,10) ## cbar for shading
        #colorbar
        levs_ticks=np.arange(0.0,121,20)
        ###Fills up the map with colors for SigmaEq
        CS1 = ax.contourf(X/1000, Y/1000, (SigmaI_Cav-SigmaI_NoCav)*1000, clevs, cmap=cmap,extend='both')
        ax.plot(xc / 1000, yc / 1000, color='k', linewidth=2)
        ###Below we remove colors that are outside of glacier contour
        clippath = mpltPath(np.c_[xc/1000, yc/1000])
        patch = PathPatch(clippath, facecolor='none')
        # ax = plt.gca()
        ax.add_patch(patch)
        for c in [CS1]:
            c.set_clip_path(patch)
        ##plot crevasses as continuous line
        for crev_num in np.arange(Df_Crevasses['Crevasse Number'].min(), Df_Crevasses['Crevasse Number'].max() + 1):
            ###get proper points
            Df_plot = Df_Crevasses[Df_Crevasses['Crevasse Number'] == crev_num]
            if Df_plot['IsCircular'].all():  ##different colors for circular crevasses and other crevasses
                col = Col_Crevasses_Circ
            else:
                col = Col_Crevasses_Other
            ax.plot(Df_plot['X'].values/1000, Df_plot['Y'].values/1000, color=col, linestyle='-', linewidth=2)
        ###Plot cavity contour
        ax.plot(cavity_contour[:, 0]/1000,cavity_contour[:, 1]/1000,color=Col_Cavity,linestyle='-',linewidth=2.5)
        ###Make a horizontal color bar common to all subplots
        # Create a colorbar axis below all subplots
        cbar_ax = fig.add_axes([0.2, 0.08, 0.6, 0.03])  # [left, bottom, width, height]
        # Create a single colorbar for all subplots
        cbar = plt.colorbar(CS1, cax=cbar_ax, ticks=levs_ticks, orientation='horizontal')
        cbar.set_label(r'$\sigma_\mathrm{I}$ [kPa]', fontsize=20)
        cbar.ax.tick_params(labelsize=18)  # Adjust tick label size


    ###Show map
    plt.show()

