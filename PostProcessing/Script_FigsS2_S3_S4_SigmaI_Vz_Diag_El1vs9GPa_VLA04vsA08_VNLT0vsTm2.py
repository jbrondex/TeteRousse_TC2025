"""
@author: jbrondex

Description:
------------
This file aims at doing three plots (one per constitutive law: elastic, linear viscous, non-linear Glen),
of the diagnostic simulation considering the empty cavity only.
Each plot is made of two lines of three subplots each. First line is surface max. principal stress and second line
is vertical displacement/velocity. First column is E= 1GPa for elastic, A=0.8 for viscous line, T=0°C for non_linear Glen.
Second column is E=9GPa for elastic, A=0.4 for viscous line, T=-2°C for non-linear Glen.
Third column is difference between first and second column.
"""

################################################################################
# import #####
################################################################################
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cmx
from matplotlib.path import Path as mpltPath
from matplotlib.patches import PathPatch
from mpl_toolkits.axes_grid1 import make_axes_locatable

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

def scientific_formatter(x, pos):
    """ Formats ticks in scientific notation with x10^exp at each tick. """
    if x == 0:
        return r"$0$"
    exp = int(np.floor(np.log10(abs(x))))  # Get exponent
    coeff = x / 10**exp  # Normalize value
    return r"${:.1f} \times 10^{{{}}}$".format(coeff, exp)

###To create colorbars for restricted number of specific subplots
def add_horizontal_colorbar(fig, CS, ax_list, label, levs_ticks, x_offset=0.0, y_offset=-0.08, sci_notation=False):
    """ Adds a single horizontal colorbar below a group of subplots. """
    cax = fig.add_axes([ax_list[0].get_position().x0 + x_offset,  # Left
                        ax_list[0].get_position().y0 + y_offset,  # Below the lowest ax in group
                        ax_list[-1].get_position().x1 - ax_list[0].get_position().x0,  # Width covering all group axes
                        0.02])  # Height of colorbar
    cbar = fig.colorbar(CS, cax=cax, orientation="horizontal", ticks=levs_ticks)
    cbar.set_label(label, fontsize=20)
    cbar.ax.tick_params(labelsize=18)
    # Apply scientific notation if requested
    if sci_notation:
        cbar.ax.xaxis.set_major_formatter(ticker.FuncFormatter(scientific_formatter))

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

    Colormap_for_dispmap = 'inferno' ###Colorm map to choose for map of stress
    cmap_disp = cmx.get_cmap(Colormap_for_dispmap )
    # cmap_disp = colors.LinearSegmentedColormap.from_list("inferno_darker", cmap_disp(np.linspace(0.2, 1, 256)))
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
    ###We first open the elastic data
    Data_Simu_El = pd.DataFrame()  ##Initialize an empty dataframe
    for E in ['1GPa', '9GPa']:
        filename_el = 'SurfaceOutput_DiagElastic_Poisson03_E{}_EmptyCavity_NoDispLat_.dat'.format(E)
        print('Opening file:', filename_el)
        Data_Simu_El_tmp = pd.read_csv(Pathroot_SimuOutput.joinpath(filename_el), names=Col_Names_Elastic, delim_whitespace=True)
        Data_Simu_El_tmp['E'] = E
        Data_Simu_El = pd.concat([Data_Simu_El, Data_Simu_El_tmp], ignore_index=True)
    ###We then open the viscous line data
    Data_Simu_VL = pd.DataFrame()
    for A in ['04','08']:
        filename_vl = 'SurfaceOutput_DiagStokes_Line_EmptyCavity_A{}_NoSliding_NoDispLat.dat'.format(A)
        print('Opening file:', filename_vl)
        Data_Simu_VL_tmp = pd.read_csv(Pathroot_SimuOutput.joinpath(filename_vl), names=Col_Names_Viscous, delim_whitespace=True)
        Data_Simu_VL_tmp['A'] = A
        Data_Simu_VL = pd.concat([Data_Simu_VL, Data_Simu_VL_tmp], ignore_index=True)
    ###We finally open the non linear Glen data
    Data_Simu_VNL = pd.DataFrame()
    for T in ['T0','Tm2']:
        filename_vnl = 'SurfaceOutput_DiagStokes_NonLine_EmptyCavity_{}_NoSliding_NoDispLat.dat'.format(T)
        print('Opening file:', filename_vnl)
        Data_Simu_VNL_tmp = pd.read_csv(Pathroot_SimuOutput.joinpath(filename_vnl), names=Col_Names_Viscous, delim_whitespace=True)
        Data_Simu_VNL_tmp['T'] = T
        Data_Simu_VNL = pd.concat([Data_Simu_VNL, Data_Simu_VNL_tmp], ignore_index=True)
    ###Create a common refined regular grid on which all fields for all cases will be interpolated
    x = np.arange(np.floor(np.min(Data_Simu_El['X'])) - 10,np.floor(np.max(Data_Simu_El['X'])) + 10, 0.1)
    y = np.arange(np.floor(np.min(Data_Simu_El['Y'])) - 10, np.floor(np.max(Data_Simu_El['Y'])) + 10, 0.1)
    X, Y = np.meshgrid(x, y)

    #######################
    ####  DO THE PLOTS:####
    #######################
    ####One plot per constitutive law (elastic, viscous line, Glen non line) but lot of common features
    for k,law in enumerate(['El','VL','VNL']):
        df_soft =pd.DataFrame()
        df_stiff=pd.DataFrame()
        ###Start by setting features that are specific to each plot
        if law == 'El':
            df_soft =  Data_Simu_El[Data_Simu_El['E']=='1GPa']
            df_soft = df_soft.sort_values(by='NodeNumber').reset_index(drop=True)
            df_stiff= Data_Simu_El[Data_Simu_El['E']=='9GPa']
            df_stiff = df_stiff.sort_values(by='NodeNumber').reset_index(drop=True)
            Deformation = 'DispZ'
            Titles = [r'$E = 1~\mathrm{GPa}$', r'$E = 9~\mathrm{GPa}$', 'Difference']
            ColBarLabel = r'$u_\mathrm{z}~\mathrm{[mm]}$'
        elif law == 'VL':
            df_soft =  Data_Simu_VL[Data_Simu_VL['A']=='08']
            df_soft = df_soft.sort_values(by='NodeNumber').reset_index(drop=True)
            df_stiff = Data_Simu_VL[Data_Simu_VL['A']=='04']
            df_stiff = df_stiff.sort_values(by='NodeNumber').reset_index(drop=True)
            Deformation = 'W'
            Titles = [r'$A = 0.8~\mathrm{MPa^{-1}~a^{-1}}$', r'$A = 0.4~\mathrm{MPa^{-1}~a^{-1}}$', 'Difference']
            ColBarLabel = r'$v_\mathrm{z}~\mathrm{[m~a^{-1}]}$'
        elif law=='VNL':
            df_soft =  Data_Simu_VNL[Data_Simu_VNL['T']=='T0']
            df_soft = df_soft.sort_values(by='NodeNumber').reset_index(drop=True)
            df_stiff = Data_Simu_VNL[Data_Simu_VNL['T']=='Tm2']
            df_stiff = df_stiff.sort_values(by='NodeNumber').reset_index(drop=True)
            Deformation = 'W'
            Titles = [r'$T = 0~\mathrm{°C}$', r'$T = -2~\mathrm{°C}$', 'Difference']
            ColBarLabel = r'$v_\mathrm{z}~\mathrm{[m~a^{-1}]}$'
        ###Prepare the subplot
        nrows = 2  ###Sigma_I, Uz
        ncols = 3  ###1GPa, 9GPa, Diff
        fig, axes = plt.subplots(nrows, ncols, figsize=(15, 12), sharex=True, sharey=True, constrained_layout=True)  #, gridspec_kw={'hspace': -0.15})
        # Store contour plots for later colorbar assignment
        CS_soft_stress = None  # Colorbar for subplots [0,0] and [0,1]
        CS_diff_stress = None  # Colorbar for subplot [0,2]
        CS_stiff_stress = None  # Colorbar for subplots [1,0], [1,1], and [1,2]
        CS_soft_disp = None  # Colorbar for subplots [0,0] and [0,1]
        CS_diff_disp = None  # Colorbar for subplot [0,2]
        CS_stiff_disp = None  # Colorbar for subplots [1,0], [1,1], and [1,2]
        for i, field in enumerate(['SigmaI', Deformation]):##loop on lines
            ###Calculate difference between soft and stiff before interpolation
            DiffField_SoftMinusStiff = df_soft[field] - df_stiff [field]
            ###Interpolate the field of interest (displacement or stress on grid)
            Field_soft_interp = []
            Field_stiff_interp = []
            DiffField_SoftMinusStiff_interp = []
            Field_soft_interp = Interpolate_field(df_soft, field, X, Y)
            Field_stiff_interp = Interpolate_field(df_stiff, field, X, Y)
            DiffField_SoftMinusStiff_interp = griddata((df_soft['X'].values, df_soft['Y'].values), DiffField_SoftMinusStiff.values, (X, Y), rescale=True)
            for j in range(ncols): ##loop on columns
                ##get proper ax
                ax = axes[i,j]
                ax.set_xlabel(r'X [km]', fontsize=22)
                ax.set_ylabel(r'Y [km]', fontsize=22)
                ax.tick_params(labelsize=18)  # fontsize of the tick labels
                ax.set_xlabel("")
                ax.set_ylabel("")
                ax.set_xticklabels([])  # Remove x labels
                ax.set_yticklabels([])  # Remove y labels
                ax.grid(True)
                ax.grid(alpha=0.5)
                ###Force x and y axis to same scale
                ax.set_aspect('equal', adjustable='box')
                if i == 0: ##titles on first line only
                    ax.set_title(Titles[j], fontsize=21, weight='bold')
                #################
                ###Do the plot###
                #################
                ####What is specific to each subplot
                if j==0: ##Softer case
                    if i==0: ## SigmaI
                        clevs = np.arange(-0.001, 121, 10)  ## cbar for shading
                        levs_ticks_stress = np.arange(0.0, 121, 20)
                        CS_soft_stress = ax.contourf(X/1000, Y/1000, Field_soft_interp*1000, clevs, cmap=cmap,extend='both')
                        C1 = CS_soft_stress
                    else: ##Deformation Z
                        if law == 'El': ##elastic case, deformation in mm
                            clevs = np.arange(-14.0001, 0.00001, 1)  ## cbar for shading
                            levs_ticks_disp = np.arange(-14, 0.0001, 2) ###Carreful : mm of displacement for this case !!!
                            CS_soft_disp = ax.contourf(X/1000, Y/1000, Field_soft_interp*1000, clevs, cmap=cmap_disp,extend='both')
                        elif law == 'VL': ###Viscous deformation in m/a
                            clevs = np.arange(-3.0, 0.0001, 0.25)  ## cbar for shading
                            levs_ticks_disp = np.arange(-3.0, 0.0001,0.5)  ###Carreful : mm of displacement for this case !!!
                            CS_soft_disp = ax.contourf(X / 1000, Y / 1000, Field_soft_interp, clevs, cmap=cmap_disp, extend='both')
                        elif law == 'VNL': ###Viscous deformation in m/a
                            clevs = np.arange(-6.0, 0.0001, 0.5)  ## cbar for shading
                            levs_ticks_disp = np.arange(-6.0, 0.0001,1.0)  ###Carreful : mm of displacement for this case !!!
                            CS_soft_disp= ax.contourf(X / 1000, Y / 1000, Field_soft_interp, clevs, cmap=cmap_disp, extend='both')
                        C1 = CS_soft_disp
                elif j ==1: ###Stiffer case
                    if i==0: ## SigmaI
                        clevs = np.arange(-0.01, 121, 10)  ## cbar for shading
                        CS_stiff_stress = ax.contourf(X/1000, Y/1000, Field_stiff_interp*1000, clevs, cmap=cmap,extend='both')
                        C1 = CS_stiff_stress
                    else:
                        if law == 'El': ##elastic case, deformation in mm
                            clevs = np.arange(-14, 0.00001, 1)  ## cbar for shading
                            levs_ticks_disp = np.arange(-14, 0.0001, 2)
                            CS_stiff_disp = ax.contourf(X/1000, Y/1000, Field_stiff_interp*1000, clevs, cmap=cmap_disp,extend='both')
                        elif law == 'VL':  ###Viscous deformation in m/a
                            clevs = np.arange(-3.0, 0.0001, 0.25)  ## cbar for shading
                            levs_ticks_disp = np.arange(-3.0, 0.0001,0.5)  ###Carreful : mm of displacement for this case !!!
                            CS_stiff_disp = ax.contourf(X / 1000, Y / 1000, Field_stiff_interp, clevs, cmap=cmap_disp,extend='both')
                        elif law == 'VNL':  ###Viscous deformation in m/a
                            clevs = np.arange(-6.0, 0.0001, 0.5)  ## cbar for shading
                            levs_ticks_disp = np.arange(-6.0, 0.0001,1.0)  ###Carreful : mm of displacement for this case !!!
                            CS_stiff_disp = ax.contourf(X / 1000, Y / 1000, Field_stiff_interp, clevs, cmap=cmap_disp,extend='both')
                        C1 = CS_stiff_disp
                else: ###Difference
                    if i==0: ## SigmaI
                        clevs = np.arange(-1e-4, 1e-4, 5e-6)  ## cbar for shading
                        levs_ticks_diff_stress = [-1e-4, 0, 1e-4]
                        CS_diff_stress = ax.contourf(X/1000, Y/1000, DiffField_SoftMinusStiff_interp*1000, clevs, cmap=cmap,extend='both')
                        C1 = CS_diff_stress
                    else:
                        if law == 'El': ##elastic case, deformation in mm
                            clevs = np.arange(-14, 0.00001, 1)
                            levs_ticks_disp = np.arange(-14, 0.0001, 2) ###Carreful : mm of displacement for this case !!!
                            CS_diff_disp = ax.contourf(X/1000, Y/1000, DiffField_SoftMinusStiff_interp*1000, clevs, cmap=cmap_disp,extend='both')
                        elif law == 'VL':  ###Viscous deformation in m/a
                            clevs = np.arange(-3.0, 0.0001, 0.25)  ## cbar for shading
                            levs_ticks_disp = np.arange(-3.0, 0.0001, 0.5)  ###Carreful : mm of displacement for this case !!!
                            CS_diff_disp = ax.contourf(X / 1000, Y / 1000, DiffField_SoftMinusStiff_interp, clevs, cmap=cmap_disp,extend='both')
                        elif law == 'VNL':  ###Viscous deformation in m/a
                            clevs = np.arange(-6.0, 0.0001, 0.5)  ## cbar for shading
                            levs_ticks_disp = np.arange(-6.0, 0.0001, 1.0)  ###Carreful : mm of displacement for this case !!!
                            CS_diff_disp = ax.contourf(X / 1000, Y / 1000, DiffField_SoftMinusStiff_interp, clevs, cmap=cmap_disp,extend='both')
                        C1 = CS_diff_disp
                ax.plot(xc / 1000, yc / 1000, color='k', linewidth=2)
                ###Below we remove colors that are outside of glacier contour
                clippath = mpltPath(np.c_[xc/1000, yc/1000])
                patch = PathPatch(clippath, facecolor='none')
                # ax = plt.gca()
                ax.add_patch(patch)
                for c in [C1]:
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
        # Add Colorbar for [0,0] and [0,1]
        add_horizontal_colorbar(fig, CS_soft_stress, [axes[0, 0], axes[0, 1]], r'$\sigma_\mathrm{I}$ [kPa]', levs_ticks_stress, x_offset=-0.03505, y_offset=-0.065)
        # Add Colorbar for [0,2] (Difference)
        add_horizontal_colorbar(fig, CS_diff_stress, [axes[0, 2]], "Difference [kPa]", levs_ticks_diff_stress, x_offset=0.05 ,y_offset=-0.065, sci_notation=True)
        # Add Colorbar for [1,0], [1,1], and [1,2]
        add_horizontal_colorbar(fig, CS_stiff_disp, [axes[1, 0], axes[1, 1], axes[1, 2]], ColBarLabel, levs_ticks_disp, x_offset=-0.01, y_offset=-0.137)

        plt.show()