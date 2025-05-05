"""
@author: jbrondex

Description:
------------
This file aims at plotting two subplots corresponding to the Maxwell time at the surface
and a vertical profile along y = y_borehole 2

"""

################################################################################
# import #####
################################################################################
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cmx
from matplotlib.path import Path as mpltPath
from matplotlib.patches import PathPatch
from matplotlib.colors import LogNorm

from pathlib import Path
from datetime import datetime, date, timedelta ###To handle dates
from matplotlib.dates import MonthLocator, DayLocator, HourLocator, DateFormatter, drange
import numpy as np
from scipy.interpolate import griddata
from scipy.interpolate import interp1d
from scipy.optimize import fsolve
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
def Interpolate_field(df, field_name, x, y): ###Returned field interpolated from mesh grid to point (x,y)
    # x, y, and z coordinates of the dataframe
    xi = df['X'].values
    yi = df['Y'].values
    field = df[field_name].values
    # Interpolate the value of field
    result = griddata((xi, yi), field, (x, y), rescale=True)
    return result

def Interpolate_field_slice(df, field_name, x, z): ###Returned field interpolated from mesh grid to point (x,z) of transect
    # x, y, and z coordinates of the dataframe
    xi = df['X'].values
    zi = df['Z'].values
    field = df[field_name].values
    # Interpolate the value of field
    result = griddata((xi, zi), field, (x, z), rescale=True)
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
    ### Path to and Name of file to be loaded for the surface
    Pathroot_SimuOutput = Path('/home/brondexj/BETTIK/TeteRousse/MyTeterousse_GeoGag/ScalarOutput/.')
    Col_Names_surf = ['Timestep', 'BC', 'NodeNumber', 'X', 'Y', 'Z', 'U', 'V', 'W', 'Pressure', 'SigmaIII','SigmaII', 'SigmaI', 'temperature', 'Sxx', 'Syy', 'Szz', 'Sxy', 'Sxz', 'Syz']
    filename_surf = 'SurfaceOutput_DiagStokes_NonLine_EmptyCavity_Tmap_NoSliding_NoDispLat.dat'
    ### Path to and Name of file to be loaded for the vertical slice
    Pathroot_SimuOutput_Slice = Path('/home/brondexj/BETTIK/TeteRousse/MyTeterousse_GeoGag/RunDiagElastic_NoDamage/ScalarOutput/.')
    Col_Names_slice = ['NodeNumber', 'SigmaIII','SigmaII', 'SigmaI', 'X', 'Y', 'Z', 'Pressure', 'Sxx', 'Syy', 'Szz', 'Sxy', 'Sxz', 'Syz', 'temperature','U', 'V', 'W']
    filename_slice = 'SurfaceOutput_DiagStokes_NonLine_EmptyCavity_Tmap_NoSliding_NoDispLat_VerticalProfil_y2105054.csv'
    ###Which Young's modulus to consider for the plot ?
    E = 9000 ### In MPa !!!!
    ### Linear or non linear Glen ?
    n = 3
    ### At which y do we want the vertical slice ? y=2105054.455 is position of borehole 2
    slice_y = 2105054.455
    ###Parameters of plot
    Col_Crevasses_Circ = '#FF0000' ###'#800080' ###Color for representation of circular crevasses
    Col_Crevasses_Other = '#800080' ###Color for representation of non circular crevasses
    Col_Cavity = '#00FFFF'###color for representation of cavity
    Col_Transect = 'k'#'#000080'###color for representation of transect

    Marker_Size = 40 ### Size of markers for crevasses on vertical slice

    Colormap_for_MaxwellTime = 'cividis_r' ###Colorm map to choose for map of stress
    cmap = cmx.get_cmap(Colormap_for_MaxwellTime  )
    cmap = colors.LinearSegmentedColormap.from_list("viridis_darker", cmap(np.linspace(0.2, 1, 256)))

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
    ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ###Load surface and bottom DEM
    Filename_zs = 'MNTSurf2011_lidar.dat'
    df_zs_DEM = pd.read_csv(Pathroot_Obs.joinpath(Filename_zs), names=['X', 'Y', 'Z'], delim_whitespace=True)
    Filename_zb = 'MNTBedSup2014.dat'
    df_zb_DEM = pd.read_csv(Pathroot_Obs.joinpath(Filename_zb), names=['X', 'Y', 'Z'], delim_whitespace=True)
    #########################################
    ####  OPEN OUTPUT OF THE SIMULATIONS:####
    #########################################
    #### For the surface directly from the Save Line solver
    Df_Simu_Surf = pd.read_csv(Pathroot_SimuOutput.joinpath(filename_surf), names=Col_Names_surf, delim_whitespace=True)
    ###Create additionnal columns to store
    ##Fluidity:
    Df_Simu_Surf['A'] = 2.42736e-02*365.25*24*60*60*1.0e18*np.exp(-115.0e3/(8.31*(273.15+Df_Simu_Surf['temperature'])))
    #Df_Simu_Surf['A'] = 2.42736e-02*365.25*24*60*60*1.0e18*np.exp(-115.0e3/(8.31*272.15))
    ###Second Invariant of the Deviatoric Stress Tensor
    Df_Simu_Surf['Is2'] = np.sqrt(0.5*((Df_Simu_Surf['Sxx'] +Df_Simu_Surf['Pressure'])**2+(Df_Simu_Surf['Syy'] +Df_Simu_Surf['Pressure'])**2+(Df_Simu_Surf['Szz'] +Df_Simu_Surf['Pressure'])**2
                                          +2*Df_Simu_Surf['Sxy']**2+2*Df_Simu_Surf['Sxz']**2+2*Df_Simu_Surf['Syz']**2))
    ###Eta as calculated in Elmer code
    Df_Simu_Surf['Eta'] = 0.5*(1/Df_Simu_Surf['A'])*Df_Simu_Surf['Is2']**(1-n)
    ###Maxwell time in year
    Df_Simu_Surf['MaxwellTime_Year'] = (2*Df_Simu_Surf['Eta']*1.3)/E

    ### Same thing for the vertical slice from the spreadsheet extracted in paraview
    Df_Simu_Slice= pd.read_csv(Pathroot_SimuOutput_Slice.joinpath(filename_slice), names=Col_Names_slice, sep=",",skiprows=1)
    ##Fluidity:
    Df_Simu_Slice['A'] = 2.42736e-02*365.25*24*60*60*1.0e18*np.exp(-115.0e3/(8.31*(273.15+Df_Simu_Slice['temperature'])))
    ###Second Invariant of the Deviatoric Stress Tensor
    Df_Simu_Slice['Is2'] = np.sqrt(0.5*((Df_Simu_Slice['Sxx'] +Df_Simu_Slice['Pressure'])**2+(Df_Simu_Slice['Syy'] +Df_Simu_Slice['Pressure'])**2+(Df_Simu_Slice['Szz'] +Df_Simu_Slice['Pressure'])**2
                                          +2*Df_Simu_Slice['Sxy']**2+2*Df_Simu_Slice['Sxz']**2+2*Df_Simu_Slice['Syz']**2))
    ###Eta as calculated in Elmer code
    Df_Simu_Slice['Eta'] = 0.5*(1/Df_Simu_Slice['A'])*Df_Simu_Slice['Is2']**(1-n)
    ###Maxwell time in year
    Df_Simu_Slice['MaxwellTime_Year'] = (2*Df_Simu_Slice['Eta']*1.3)/E

    #########################################################################
    ###     PLOT MAP OF MAXWELL TIME AT SURFACE AND OVER VERTICAL PROFILE ###
    #########################################################################
    ###~~~~~~~Prepare the figure 1 : subplot with the surface eq stress for the four criterion tested
    ###Prepare the subplot
    nrows = 1  ###
    ncols = 2  ###Surf, Vertical Profile
    fig1, axes = plt.subplots(nrows, ncols, figsize=(15, 12), constrained_layout=True)  # , gridspec_kw={'hspace': -0.15})
    # plt.subplots_adjust(bottom=0.2,wspace=0.1)
    ####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####
    ####        START WITH SUBPLOT AT THE SURFACE           ####
    ####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####
    ax = axes[0]  # Select the axis in the grid
    ax.set_xlabel(r'X [km]', fontsize=18)
    ax.set_ylabel(r'Y [km]', fontsize=18)
    ax.tick_params(labelsize=16)  # fontsize of the tick labels
    ax.grid(True)
    ax.grid(alpha=0.5)
    ###Force x and y axis to same scale
    ax.set_aspect('equal', adjustable='box')
    # ax.set_title(Title, fontsize=21, weight='bold')

    ###Create refined regular grid and interpolate field over grid
    xsurf = np.arange(np.floor(np.min(Df_Simu_Surf['X']))-10,np.floor(np.max(Df_Simu_Surf['X']))+10,0.1)
    ysurf = np.arange(np.floor(np.min(Df_Simu_Surf['Y']))-10,np.floor(np.max(Df_Simu_Surf['Y']))+10,0.1)
    Xsurf, Ysurf = np.meshgrid(xsurf, ysurf)
    MaxwellTime_surf = Interpolate_field(Df_Simu_Surf,'MaxwellTime_Year',Xsurf,Ysurf)
    #shading
    clevs=np.logspace(0, 4, num=9)  ## cbar for shading
    cmap=Colormap_for_MaxwellTime
    #colorbar
    levs_ticks = [1, 10, 100, 1000, 10000]
    ###Fills up the map with colors for SigmaEq
    CS1 = ax.contourf(Xsurf/1000, Ysurf/1000, MaxwellTime_surf*365.25*24, levels=clevs, cmap=cmap, norm=LogNorm(vmin=1, vmax=1e4), extend='both')
    ax.plot(xc / 1000, yc / 1000, color='k', linewidth=2)
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
            col = Col_Crevasses_Other
        ax.plot(Df_plot['X'].values/1000, Df_plot['Y'].values/1000, color=col, linestyle='-', linewidth=2)
    ###plot transect corresponding to vertical slice
    ax.axhline(y=slice_y / 1000, color=Col_Transect, linestyle='--', linewidth=2.2)
    ### Set ylim
    ax.set_ylim([
        (np.floor(np.min(Df_Simu_Surf['Y'])) - 10) / 1000,
        (np.floor(np.max(Df_Simu_Surf['Y'])) + 10) / 1000
    ])

    ####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####
    ####      NOW THE SUBPLOT FOR THE VERTICAL SLICE        ####
    ####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####
    ax = axes[1]  # Select the axis in the grid
    ax.set_xlabel(r'X [km]', fontsize=18)
    ax.set_ylabel(r'Z [m]', fontsize=18)
    ax.tick_params(labelsize=16)  # fontsize of the tick labels
    ax.grid(True)
    ax.grid(alpha=0.5)
    ###Force x and y axis to same scale
    ax.set_aspect(1/1000)  # Ensures equal unit length on both axes
    ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##########
    ###Extract a 2D vertical slice from the surf and bottom DEM -> For visualization only ###
    ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##########
    # Define the vertical slice line (a profile along y=y of borehole 2)
    tolerance_zs = 3  # Chose tolerance so that 2 points are kept above and below the slice (in y_direction) for each x value
    tolerance_zb = 2
    # Keep all DEM points close enough to the vertical slice
    slice_points_zs = df_zs_DEM[(df_zs_DEM['Y'] > slice_y - tolerance_zs) & (df_zs_DEM['Y'] < slice_y + tolerance_zs)]
    slice_points_zb = df_zb_DEM[(df_zb_DEM['Y'] > slice_y - tolerance_zb) & (df_zb_DEM['Y'] < slice_y + tolerance_zb)]
    ###Check that points that have been kept contains no more than two y values
    y_zs_DEM_lines = slice_points_zs['Y'].drop_duplicates()
    if len(y_zs_DEM_lines) !=2:
        raise ValueError("Vertical slice for zs DEM is built from more (or less) than 2 points in y-direction. Tolerance must be changed")
    else:
        y1_zs = y_zs_DEM_lines.iloc[0]
        y2_zs = y_zs_DEM_lines.iloc[1]
    ###Same thing for zb
    y_zb_DEM_lines = slice_points_zb['Y'].drop_duplicates()
    if len(y_zb_DEM_lines) !=2:
        raise ValueError("Vertical slice for zb DEM is built from more (or less) than 2 points in y-direction. Tolerance must be changed")
    else:
        y1_zb = y_zb_DEM_lines.iloc[0]
        y2_zb = y_zb_DEM_lines.iloc[1]
    ###Now build df_zs_inter where y=y_slice, x is extracted from slice_points_zs and z is interpolated linearly at y=slice_y from z at y1 and z at y2
    df_zs_interp=pd.DataFrame(columns=['X', 'Y', 'Z'])
    for index, row1 in slice_points_zs[slice_points_zs['Y']==y1_zs].iterrows():
        ###get corresponding row for y = y2
        row2 = slice_points_zs[(slice_points_zs['Y']==y2_zs) & (slice_points_zs['X']==row1['X'])].iloc[0]
        ### Check that x values of points are the same along y=y1 and along y=y2
        if not row1['X'] == row2['X']:
            raise ValueError("x values are not strictly the same on both sides of vertical slice for zs")
        ### for the current x, interpolate z at y=slice_y
        slice_zs = np.interp(slice_y, [y1_zs, y2_zs], [row1['Z'], row2['Z']])
        ### Fill up interpolated data frame
        df_zs_interp.loc[len(df_zs_interp)] = [row1['X'], slice_y, slice_zs]
    ###Same thing for zb
    df_zb_interp=pd.DataFrame(columns=['X', 'Y', 'Z'])
    for index, row1 in slice_points_zb[slice_points_zb['Y']==y1_zb].iterrows():
        ###get corresponding row for y = y2
        row2 = slice_points_zb[(slice_points_zb['Y']==y2_zb) & (slice_points_zb['X']==row1['X'])].iloc[0]
        ### Check that x values of points are the same along y=y1 and along y=y2
        if not row1['X'] == row2['X']:
            raise ValueError("x values are not strictly the same on both sides of vertical slice for zb")
        ### for the current x, interpolate z at y=slice_y
        slice_zb = np.interp(slice_y, [y1_zb, y2_zb], [row1['Z'], row2['Z']])
        ### Fill up interpolated data frame
        df_zb_interp.loc[len(df_zb_interp)] = [row1['X'], slice_y, slice_zb]
    ####Remove the parts below and above intersections
    # Interpolation functions for both lines
    f1 = interp1d(df_zb_interp['X'], df_zb_interp['Z'], kind='linear', fill_value="extrapolate")
    f2 = interp1d(df_zs_interp['X'], df_zs_interp['Z'], kind='linear', fill_value="extrapolate")
    # Define function to find intersection points (f1(x) - f2(x) = 0)
    def find_intersection(x):
        return f1(x) - f2(x)
    # Estimate intersections (adjust range if necessary)
    x_intersections = fsolve(find_intersection, [np.min(Df_Simu_Surf['X']), 948272])  # Initial guesses
    # Extract values
    x_intersection1, x_intersection2 = sorted(x_intersections)
    # Keep only the parts of the lines between the two intersections
    df_zs_interp_filtered= df_zs_interp[(df_zs_interp['X'] >= x_intersection1) & (df_zs_interp['X'] <= x_intersection2)]
    df_zb_interp_filtered = df_zb_interp[(df_zb_interp['X'] >= x_intersection1) & (df_zb_interp['X'] <= x_intersection2)]
    ### build contour glacier from zs and zb
    df_zs_interp_filtered= df_zs_interp_filtered.iloc[::-1].reset_index(drop=True)
    df_contour_slice = pd.concat([df_zb_interp_filtered,df_zs_interp_filtered])
    ###close contour
    last_row = df_contour_slice.iloc[[-1]]
    df_contour_slice = pd.concat([last_row,df_contour_slice], ignore_index=True)
    ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ### Create interpolation grid for the vertical slice and interpolate
    xslice = xsurf
    zmin = df_zb_DEM['Z'].min()
    zmax = 3270
    zslice = np.arange(zmin, zmax, 0.1)
    Xslice, Zslice = np.meshgrid(xslice, zslice)
    ###INTERPOLATE Maxwell time on grid
    MaxwellTime_slice = Interpolate_field_slice(Df_Simu_Slice,'MaxwellTime_Year',Xslice,Zslice)
    ###Fills up the map with colors for SigmaEq
    CS2 = ax.contourf(Xslice/1000, Zslice, MaxwellTime_slice*365.25*24, levels=clevs, cmap=cmap, norm=LogNorm(vmin=1, vmax=1e4), extend='both')
    ax.plot(df_contour_slice['X'].values/1000, df_contour_slice['Z'].values, color='k', linewidth=2)
    ###Below we remove colors that are outside of glacier contour
    clippath = mpltPath(np.c_[df_contour_slice['X'].values/1000, df_contour_slice['Z'].values])
    patch = PathPatch(clippath, facecolor='none')
    ax.add_patch(patch)
    for c in [CS2]:
        c.set_clip_path(patch)
    ### Add marks for crevasses crossed by transect
    ###circular crevasses (coordinates are found manually)
    ax.scatter(947953.5 / 1000, 3162.5, color=Col_Crevasses_Circ, marker='o', s=Marker_Size, zorder=3)
    ax.scatter(948091.4 / 1000, 3193, color=Col_Crevasses_Circ, marker='o', s=Marker_Size, zorder=3)
    ax.scatter(948103.5 / 1000, 3196.8, color=Col_Crevasses_Circ, marker='o', s=Marker_Size, zorder=3)
    ### Other crevasse
    ax.scatter(948177.5 / 1000, 3224, color=Col_Crevasses_Other, marker='o', s=Marker_Size, zorder=3)
    ### The vertical extent of the subplot has to be adjusted so that both subplots have same height while keeping 1/1 aspect ratio between axes
    PlotWidth_AlongY = (np.floor(np.max(Df_Simu_Surf['Y']))+10) - (np.floor(np.min(Df_Simu_Surf['Y']))-10)
    PlotHeight_AlongZ = zmax - zmin
    Diff_Width_Height = PlotWidth_AlongY - PlotHeight_AlongZ
    ax.set_ylim([zmin-Diff_Width_Height/2, zmax + Diff_Width_Height/2])
    ####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####
    ####      COMMON TO THE TWO SUBPLOTS        ####
    ####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####
    # Add a common colorbar to all subplots below the figure
    cbar_ax = fig1.add_axes([0.15, 0.082, 0.7, 0.027])  # [left, bottom, width, height]
    fig1.colorbar(CS1, ticks=levs_ticks, cax=cbar_ax, orientation='horizontal', label=r'Maxwell time [h]')
    cbar_ax.set_xticklabels([r"$10^0$", r"$10^1$", r"$10^2$", r"$10^3$", r"$10^4$"])
    ###Show map
    plt.show()




