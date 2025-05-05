"""
@author: jbrondex

Description:
------------
This file aims at building a 3D Temperature field of TeteRousse from borehole measurement performed the 11 September 2010
First we interpolate T field from borehole measurements, and where required we extrapolate using insights from
of Gilbert 2012
"""

################################################################################
# import #####
################################################################################
import matplotlib
from matplotlib.lines import lineStyles

matplotlib.use('TkAgg',force=True)
from matplotlib import pyplot as plt
# import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.cm as cmx
from fontTools.unicodedata import block
from matplotlib.path import Path as mpltPath
from matplotlib.patches import PathPatch

from pathlib import Path
from datetime import datetime, date, timedelta ###To handle dates
from matplotlib.dates import MonthLocator, DayLocator, HourLocator, DateFormatter, drange
import numpy as np
from scipy.spatial import cKDTree
from scipy.interpolate import interp1d
from scipy.interpolate import griddata
from scipy.optimize import fsolve
from matplotlib.ticker import FormatStrFormatter
from matplotlib.dates import date2num
import matplotlib.gridspec as gridspec

import pandas as pd  ###To treat the csv files

from matplotlib import ticker

##Function below is used to sort point of contour clockwise
def sort_clockwise(list_of_xy_coords):
    cx, cy = list_of_xy_coords.mean(0)
    x, y = list_of_xy_coords.T
    angles = np.arctan2(x-cx, y-cy)
    indices = np.argsort(angles)
    return list_of_xy_coords[indices]

def Interpolate_field(df, field_name, x, z): ###Returned field interpolated from mesh grid to point (x,y)
    # x, y, and z coordinates of the dataframe
    xi = df['X'].values
    zi = df['Sensor_Altitude'].values
    field = df[field_name].values
    # Interpolate the value of field
    result = griddata((xi, zi), field, (x, z), rescale=True)
    return result

# Function to truncate colormap
def truncate_colormap(cmap, min_val=0.0, max_val=1.0, n_colors=256):
    """Create a truncated colormap from an existing one."""
    new_cmap = mcolors.LinearSegmentedColormap.from_list(
        'truncated_cmap', cmap(np.linspace(min_val, max_val, n_colors))
    )
    return new_cmap

###Function to replace NaN of a griddata field obtained by interpolation by value of closest neighbour
def fill_nan_with_nearest(T, X, Z):
    # Get valid (non-NaN) points
    valid_mask = ~np.isnan(T)
    valid_points = np.column_stack([X[valid_mask], Z[valid_mask]])
    valid_values = T[valid_mask]

    # Get NaN points
    nan_mask = np.isnan(T)
    nan_points = np.column_stack([X[nan_mask], Z[nan_mask]])

    # Build KDTree and query nearest neighbor
    tree = cKDTree(valid_points)
    _, nearest_idx = tree.query(nan_points)

    # Replace NaNs with nearest neighbor values
    T[nan_mask] = valid_values[nearest_idx]

    return T

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
    #### Parameters of boreholes for the 11 September 2010
    Borehole_List = [2, 4, 5, 10, 13, 17, 18]
    Borehole_Number_ForVerticalSlice = 2 ###Number of Borehole used as reference for vertical slice
    #### Do we want to clip the extrapolated T map ? If yes then set IsClip to True
    IsClip = True
    #### Parameter for plot of T map
    Tmin=-3.0
    Tmax= 0.01
    lev_ticks_range = 0.5  ##one tick every ?
    color_resolution = 10  ##number of color shades per level ticks
    ###List of color codes from colder to warmer used in Fig. 5 of Gilbert 2012
    Colors_Gilbert2012 = ['#002fdd','#0048ff','#008cff','#00cfff','#00fff5','#00fead','#98fc54','#fffa00','#ffb700','#ff6b00','#ff0000', '#d30000' ]
    ##################################################################
    ####  OPEN DATA CORRESPONDING TO GEOMETRY OF THE GLACIER     ####
    #################################################################
    ###Load glacier contour
    Pathroot_Obs = Path('/home/brondexj/BETTIK/TeteRousse/MyTeterousse_GeoGag/Data/.')
    Filename_GlacierContour = 'Contour_TeteRousse2012.dat'
    Df_GlacierContour = pd.read_csv(Pathroot_Obs.joinpath(Filename_GlacierContour), names=['X', 'Y'], delim_whitespace=True)
    xc = Df_GlacierContour['X'].values
    xc = np.append(xc, xc[0])  ##step required to close the contour
    yc = Df_GlacierContour['Y'].values
    yc = np.append(yc, yc[0])  ##step required to close the contour
    ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ###Open Data corresponding to grounded mask
    Pathroot_GM = Path('/home/brondexj/BETTIK/TeteRousse/MyTeterousse_GeoGag/ScalarOutput/.')
    Filename_GM = 'GroundedLine_CLEANED.dat'
    Col_Names_GM = ['Timestep', 'BC', 'NodeNumber', 'X', 'Y', 'Z', 'GM']
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

    ###############################################
    ####   NOW START THE REQUIRED PROCESSING   ####
    ###############################################
    ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ###OPEN DATA CORRESPONDING TO BOREHOLES
    Pathroot_T = Path('/home/brondexj/BETTIK/TeteRousse/MyTeterousse_GeoGag/PostProcessing/Raw_Temperature_TeteRousse/.')
    Filename_Info_T = 'TemperatureTR_Info.csv'
    df_Info_T = pd.read_csv(Pathroot_T.joinpath(Filename_Info_T), sep=";")
    ### do the loop on boreholes and create a single dataframe with col [Borehole, x, y, zsurf, profondeur, T]
    Col_Names_Borehole = ['Depth', 'T']
    df_AllBoreholes = pd.DataFrame()
    for i, Borehole in enumerate(Borehole_List):
        x_Borehole = df_Info_T[df_Info_T['drilling_id']==Borehole]['x_lat'].values[0]
        y_Borehole = df_Info_T[df_Info_T['drilling_id']==Borehole]['y_lon'].values[0]
        zs_Borehole = df_Info_T[df_Info_T['drilling_id']==Borehole]['zs_m'].values[0]
        Filename_Borehole = '2010-09-11/TemperatureTR_{}_2010-09-11.csv'.format(str(Borehole))
        df_SingleBorehole = pd.DataFrame()
        df_SingleBorehole = pd.read_csv(Pathroot_T.joinpath(Filename_Borehole), names=Col_Names_Borehole, sep=";", skiprows=1)
        df_SingleBorehole.insert(0, 'Borehole_ID', Borehole)
        df_SingleBorehole.insert(1, 'X', x_Borehole)
        df_SingleBorehole.insert(2, 'Y', y_Borehole)
        df_SingleBorehole.insert(3, 'Zs', zs_Borehole)
        ###  ###Altitude of each thermistance is required:
        df_SingleBorehole['Sensor_Altitude']= zs_Borehole+df_SingleBorehole['Depth']
        df_AllBoreholes = pd.concat([df_AllBoreholes,df_SingleBorehole], ignore_index=True)
    ###Remove line with NaN (there is a sensor missing in borehole 5)
    df_AllBoreholes.dropna(inplace=True)
    ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ### Build a temperature map on a vertical slice by interpolating measurements
    ###Create refined regular grid and interpolate field over grid
    xmin = xc.min()
    xmax = xc.max()
    zmin = df_zb_DEM['Z'].min()
    zmax = 3270
    x = np.arange(xmin, xmax, 2)
    z = np.arange(zmin, zmax, 2)
    X, Z = np.meshgrid(x, z)
    ###INTERPOLATE T ON GRID FROM MEASUREMENTS ONLY
    Interpolated_T = Interpolate_field(df_AllBoreholes, 'T', X, Z)
    ### Create a NaN mask (1 where NaN, 0 where valid data) (to show where data are interpolated and where they are extrapolated
    nan_mask = np.isnan(Interpolated_T).astype(float)

    ##################################################################
    ###   BUILD AN EXTRAPOLATED T MAP FROM THE INTERPOLATED T MAP  ###
    ##################################################################
    ###Copy the dataframe of Boreholes
    df_AllBoreholes_Extrapolated = df_AllBoreholes.copy()
    ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ### First, correct second shallowest sensor of borehole 18 that give a strange result
    df_borehole_18 = df_AllBoreholes_Extrapolated[df_AllBoreholes_Extrapolated['Borehole_ID']==18]

    # Find the second highest Depth row as before
    df_borehole_18_sorted = df_borehole_18.sort_values(by='Depth', ascending=False).reset_index(drop=True)
    second_highest_row_index = df_borehole_18_sorted.iloc[1].name
    depth_sensor = df_borehole_18_sorted['Depth'].loc[second_highest_row_index]
    # Get the previous and next rows in the sorted dataframe
    prev_row = df_borehole_18_sorted.iloc[second_highest_row_index - 1] if second_highest_row_index - 1 >= 0 else None
    next_row = df_borehole_18_sorted.iloc[second_highest_row_index + 1] if second_highest_row_index + 1 < len(df_borehole_18_sorted) else None
    # If both previous and next rows exist, calculate the mean of the 'T' values
    if prev_row is not None and next_row is not None:
        mean_T = (prev_row['T'] + next_row['T']) / 2
        # Update the 'T' value of the second highest row
        df_AllBoreholes_Extrapolated.loc[(df_AllBoreholes_Extrapolated['Borehole_ID']==18) & (df_AllBoreholes_Extrapolated['Depth']==depth_sensor), 'T'] = mean_T
    ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ### Second, add a virtual borhehole in the upper part of the glacier (e.g 70 m higher than upper borehole) where we know the base is temperate
    x_virtualborehole1 = df_AllBoreholes_Extrapolated['X'].max() + 70
    y_virtualborehole1 = df_AllBoreholes_Extrapolated[df_AllBoreholes_Extrapolated['X'] == df_AllBoreholes_Extrapolated['X'].max()]['Y'].iloc[0]
    zs_virtualborehole1 = 3250
    zb_virtualborehole1 = 3222
    Sensor_Altitude_virtualborehole1 = np.arange(zb_virtualborehole1, zs_virtualborehole1, 5)
    Depth_virtualborehole1 = - (zs_virtualborehole1 - Sensor_Altitude_virtualborehole1)
    Virtual_Borehole1 = pd.DataFrame({'Borehole_ID': 100 * np.ones(len(Sensor_Altitude_virtualborehole1)),
                                      'X': x_virtualborehole1 * np.ones(len(Sensor_Altitude_virtualborehole1)),
                                      'Y': y_virtualborehole1 * np.ones(len(Sensor_Altitude_virtualborehole1)),
                                      'Zs': zs_virtualborehole1 * np.ones(len(Sensor_Altitude_virtualborehole1)),
                                      'Depth': Depth_virtualborehole1,
                                      'T': 0 * np.ones(len(Sensor_Altitude_virtualborehole1)),
                                      'Sensor_Altitude': Sensor_Altitude_virtualborehole1})
    ##set T to -0.25 in vicinity of the surface for virtual borehole (Gilbert 2012)
    Virtual_Borehole1.loc[Virtual_Borehole1['Depth'] > -8, 'T'] = -0.25
    ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ### Similarly, add a virtual borhehole in the lower part of the glacier (e.g 35 m lower than lower borehole) where we know ice is cold
    x_virtualborehole2 = df_AllBoreholes_Extrapolated['X'].min() - 35
    y_virtualborehole2 = df_AllBoreholes_Extrapolated[df_AllBoreholes_Extrapolated['X'] == df_AllBoreholes_Extrapolated['X'].min()]['Y'].iloc[0]
    zs_virtualborehole2 = 3145
    zb_virtualborehole2 = 3118
    Sensor_Altitude_virtualborehole2 = np.arange(zb_virtualborehole2, zs_virtualborehole2, 5)
    Depth_virtualborehole2 = - (zs_virtualborehole2 - Sensor_Altitude_virtualborehole2)
    Virtual_Borehole2 = pd.DataFrame({'Borehole_ID': 101 * np.ones(len(Sensor_Altitude_virtualborehole2)),
                                      'X': x_virtualborehole2 * np.ones(len(Sensor_Altitude_virtualborehole2)),
                                      'Y': y_virtualborehole2 * np.ones(len(Sensor_Altitude_virtualborehole2)),
                                      'Zs': zs_virtualborehole2 * np.ones(len(Sensor_Altitude_virtualborehole2)),
                                      'Depth': Depth_virtualborehole2,
                                      'T': -2.75 * np.ones(len(Sensor_Altitude_virtualborehole2)),
                                      'Sensor_Altitude': Sensor_Altitude_virtualborehole2})
    ###Add virtual borehole to true ones
    df_AllBoreholes_Extrapolated= pd.concat([Virtual_Borehole2, df_AllBoreholes_Extrapolated,Virtual_Borehole1], ignore_index=True)

    ####Then, Add T at surface of each borehole
    for Borehole in df_AllBoreholes_Extrapolated['Borehole_ID'].drop_duplicates():
        if Borehole == 100:##don't do the operation for virtual borehole
            continue
        df_Borehole = df_AllBoreholes_Extrapolated[df_AllBoreholes_Extrapolated['Borehole_ID']==Borehole]
        ### Tsurf is deduced from T closest to surface to which we remove 0.4°C to
        Tsurf = df_Borehole.loc[df_Borehole['Depth'].idxmax(),'T'] - 0.4
        new_row = pd.DataFrame({'Borehole_ID': [Borehole], 'X': [df_Borehole['X'].iloc[0]], 'Y': [df_Borehole['Y'].iloc[0]],
                             'Zs': [df_Borehole['Zs'].iloc[0]], 'Depth': [0], 'T': [Tsurf], 'Sensor_Altitude': [df_Borehole['Zs'].iloc[0]+4]})
        # Find the index of the last occurrence of the current Borehole
        last_index = df_AllBoreholes_Extrapolated[df_AllBoreholes_Extrapolated['Borehole_ID'] == Borehole].index[-1]
        # Split the DataFrame and insert the new row after last occurrence
        df_AllBoreholes_Extrapolated = pd.concat([df_AllBoreholes_Extrapolated.iloc[:last_index+1], new_row, df_AllBoreholes_Extrapolated.iloc[last_index+1:]], ignore_index=True)

    ### Replace all values above -0.25°C by 0°C to have temperate ice around cavity
    df_AllBoreholes_Extrapolated.loc[df_AllBoreholes_Extrapolated['T'] > -0.25, 'T'] = 0
    ###INTERPOLATE T ON GRID FROM MEASUREMENTS + VIRTUAL BOREHOLES + VIRTUAL SURFACE T
    Extrapolated_T = Interpolate_field(df_AllBoreholes_Extrapolated, 'T', X, Z)
    ###Above the upper virtual borehole we force T to 0°C
    # Find NaN values where X > xmax and Z < zmin
    mask1 = (X > x_virtualborehole1) & (Z < zs_virtualborehole1+20) & np.isnan(Extrapolated_T)
    # Replace NaNs with 0 in those locations
    Extrapolated_T[mask1] = 0
    ###Similarly, below the lower virtual borehole we force T to -2.6°C
    # Find NaN values where X < xmin and Z < zmin
    mask2 = (X < x_virtualborehole2) & (Z < zs_virtualborehole2) & np.isnan(Extrapolated_T)
    # Replace NaNs with -2.6 in those locations
    Extrapolated_T[mask2] = -2.6

    # ####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####
    # ####  Produce final file and save it for interpolation through Elmer Solver  #####
    # ####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####
    # #### Replace all remaining Nan by -9999 for Elmer Solver
    # Extrapolated_T = np.nan_to_num(Extrapolated_T, nan=-9999)
    # # Flatten the grids to 1D arrays
    # X_flat = X.ravel()  # Convert 2D grid to 1D array
    # Z_flat = Z.ravel()
    # T_flat = Extrapolated_T.ravel()
    # # Stack the arrays column-wise
    # T_data = np.column_stack((X_flat, Z_flat, T_flat))
    # # Save to .dat file
    # Filename_Tmap= 'VerticallyExtrapolated_2D_Temperature.dat'
    # Pathroot_Tdata = Path('/home/brondexj/BETTIK/TeteRousse/MyTeterousse_GeoGag/Data/.')
    # np.savetxt(Pathroot_Tdata.joinpath(Filename_Tmap), T_data, fmt="%.6f", delimiter=" ", comments="")

    ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##########
    ###Extract a 2D vertical slice from the surf and bottom DEM -> For visualization only ###
    ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##########
    # Define the vertical slice line (a profile along y=y of borehole 2)
    slice_y = df_AllBoreholes[df_AllBoreholes['Borehole_ID']==Borehole_Number_ForVerticalSlice]['Y'].iloc[0]  ###Vertical slice at y of borehole 2
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
    x_intersections = fsolve(find_intersection, [xmin, 948272])  # Initial guesses
    # Extract values
    x_intersection1, x_intersection2 = sorted(x_intersections)
    # Keep only the parts of the lines between the two intersections
    df_zs_interp_filtered= df_zs_interp[(df_zs_interp['X'] >= x_intersection1) & (df_zs_interp['X'] <= x_intersection2)]
    df_zb_interp_filtered = df_zb_interp[(df_zb_interp['X'] >= x_intersection1) & (df_zb_interp['X'] <= x_intersection2)]
    ### build contour glacier from zs and zb
    df_zs_interp_filtered= df_zs_interp_filtered.iloc[::-1].reset_index(drop=True)
    df_contour_glacier = pd.concat([df_zb_interp_filtered,df_zs_interp_filtered])
    ###close contour
    last_row = df_contour_glacier.iloc[[-1]]
    df_contour_glacier = pd.concat([last_row,df_contour_glacier], ignore_index=True)
    ######################################
    ####         DO THE PLOTS          ###
    ######################################
    ###PLOT MAP OF TETEROUSSE FROM ABOVE
    ###Fig 1 is the map of Teterousse, pimp style
    fig1 = plt.figure(1, figsize=(40, 20))
    plt.xlabel(r'X [km]', fontsize=21)
    plt.ylabel(r'Y [km]', fontsize=21)
    plt.tick_params(labelsize=20)  # fontsize of the tick labels
    plt.grid(True)
    plt.grid(alpha=0.5)
    # plt.xlim([947950, 948050])
    # plt.ylim([2105000,2105120])
    ###Force x and y axis to same scale
    ax = plt.gca()
    plt.plot(xc / 1000, yc / 1000, color='k', linestyle='-', linewidth=1.5)
    plt.plot(cavity_contour[:, 0] / 1000, cavity_contour[:, 1] / 1000, color='grey', linestyle='-', linewidth=1.2)
    plt.axhline(y=slice_y/ 1000, color='r', linestyle='-', linewidth=1.8)
    ax.set_aspect('equal',
                  adjustable='box')  # plt.plot(xy_sorted[:, 0]/1000,xy_sorted[:, 1]/1000,color='dimgrey',linestyle='-',linewidth=1.5)
    plt.scatter(df_AllBoreholes['X'].values / 1000, df_AllBoreholes['Y'].values / 1000, color='c', edgecolor='black',
                marker='*', s=100, linewidths=0.6)
    # Add labels to each scatter point
    for i, Borehole in enumerate(
            Borehole_List):  # Iterate through each point and label    plt.plot(xy_sorted[:, 0],xy_sorted[:, 1],color='k',linestyle='-',linewidth=1)
        plt.text(
            df_AllBoreholes[df_AllBoreholes['Borehole_ID'] == Borehole]['X'].iloc[0] / 1000,
            df_AllBoreholes[df_AllBoreholes['Borehole_ID'] == Borehole]['Y'].iloc[0] / 1000,
            Borehole,
            fontsize=9,  # Adjust label font size
            ha='right',  # Horizontal alignment (aligns text to the right of the point)
            va='bottom'
        )
    ###Fig 2: VERTICAL SLICE WITH INTERPOLATED T MAP
    fig2 = plt.figure(2, figsize=(40, 20))
    plt.xlabel(r'X [km]', fontsize=21)
    plt.ylabel(r'Z [m]', fontsize=21)
    plt.tick_params(labelsize=20)  # fontsize of the tick labels
    plt.grid(True)
    plt.grid(alpha=0.5)
    ###Set xaxis and y axis to same scale, but knowing that x is in km
    ax = plt.gca()
    ax.set_aspect(1/1000)  # Ensures equal unit length on both axes
    ###Fills up the map with colors for SigmaEq
    ##contour of T every 0.25°C
    clevs_contour = np.arange(Tmin, Tmax, lev_ticks_range/2)
    clevs_contour = np.sort(clevs_contour)  # Ensure proper ordering
    CS = plt.contour(X/ 1000, Z, Interpolated_T,clevs_contour, colors='black', linestyles='-', linewidths=0.6) ## lines (f(x) contour)
    plt.clabel(CS, CS.levels[::2], inline=True, fontsize=8) ## numbers in contours
    ###Fills up the map with colors for T
    CS1 = plt.contourf(X / 1000, Z, Interpolated_T, clevs_contour, colors=Colors_Gilbert2012)
    ###Show colorbar
    levs_ticks = np.arange(Tmin, Tmax, lev_ticks_range)
    cbar = plt.colorbar(CS1, ticks=levs_ticks, orientation='vertical', label=r'T [°C]')
    ###plot thermistance positions for each borehole
    for i, Borehole in enumerate(Borehole_List):
        plt.scatter(df_AllBoreholes[df_AllBoreholes['Borehole_ID'] == Borehole]['X'].values/1000,
                    df_AllBoreholes[df_AllBoreholes['Borehole_ID'] == Borehole]['Zs'].values+df_AllBoreholes[df_AllBoreholes['Borehole_ID'] == Borehole]['Depth'],
                    color='lime', edgecolor='black', marker='.', s=60, zorder=10)
    ###plot contour of chosen slice
    plt.plot(df_contour_glacier['X'].values/1000, df_contour_glacier['Z'].values, color='k', linestyle='-', linewidth=1.7)

    ###Fig 3: VERTICAL SLICE WITH EXTRAPOLATED T MAP
    fig3 = plt.figure(3, figsize=(40, 20))
    plt.xlabel(r'X [km]', fontsize=21)
    plt.ylabel(r'Z [m]', fontsize=21)
    plt.tick_params(labelsize=20)  # fontsize of the tick labels
    plt.grid(True)
    plt.grid(alpha=0.5)
    ###Set xaxis and y axis to same scale, but knowing that x is in km
    ax = plt.gca()
    ax.set_aspect(1/1000)  # Ensures equal unit length on both axes
    ###plot contour of chosen slice
    plt.plot(df_contour_glacier['X'].values/1000, df_contour_glacier['Z'].values, color='k', linestyle='-', linewidth=1.7)
    # plt.plot(df_zs_interp_filtered['X'].values/1000, df_zs_interp_filtered['Z'].values-10, color='k', linestyle='--', linewidth=2)
    ##contour of T every 0.25°C
    clevs_contour = np.arange(Tmin, Tmax, lev_ticks_range / 2)
    clevs_contour = np.sort(clevs_contour)  # Ensure proper ordering
    ###Fills up the map with colors for T
    CS1 = plt.contourf(X / 1000, Z, Extrapolated_T, clevs_contour, colors=Colors_Gilbert2012, extend = 'min')
    ###Plot contour lines
    CS = plt.contour(X / 1000, Z, Extrapolated_T, clevs_contour, colors='black', linestyles='-',
                     linewidths=0.6)  ## lines (f(x) contour)
    # # Select label positions manually (problem with last version of matplotlib)
  ##  label_positions = [(947.881, 3141.5), (947.9055, 3132), (947.9466, 3134.0), (948.0326, 3162), (948.1157, 3177), (948.0773, 3147.5)]  # Adjust these positions as needed
    plt.clabel(CS, CS.levels[::2], inline=True, fontsize=8) #, manual=valid_positions)  ## numbers in contours

    ###Show colorbar
    levs_ticks = np.arange(Tmin, Tmax, lev_ticks_range)
    cbar = plt.colorbar(CS1, ticks=levs_ticks, orientation='vertical', label=r'T [°C]')
    ###Plot contour between NaN and no NaN to show where data are interpolated
    CS2 = plt.contour(X / 1000, Z, nan_mask, levels=[0.5], colors='whitesmoke', linestyles='-', linewidths=2)
    ###Remove what is out of glacier contour
    ###Below we remove colors that are outside of glacier contour
    if IsClip:
        clippath = mpltPath(np.c_[df_contour_glacier['X'].values / 1000, df_contour_glacier['Z'].values])
        patch = PathPatch(clippath, facecolor='none')
        ax = plt.gca()
        ax.add_patch(patch)
        for c in [CS]:
            c.set_clip_path(patch)
        for c in [CS1]:
            c.set_clip_path(patch)
        for c in [CS2]:
            c.set_clip_path(patch)
    ###plot thermistance positions for each borehole (True Borehole)
    for i, Borehole in enumerate(Borehole_List):
        plt.scatter(df_AllBoreholes[df_AllBoreholes['Borehole_ID'] == Borehole]['X'].values / 1000,
                    df_AllBoreholes[df_AllBoreholes['Borehole_ID'] == Borehole]['Zs'].values +
                    df_AllBoreholes[df_AllBoreholes['Borehole_ID'] == Borehole]['Depth'],
                    color='lime', edgecolor='black', marker='.', s=60, zorder=10)
    ### If you want to show the interpolation grid, uncomment following:
    # plt.scatter(X_flat/1000,Z_flat,color='k', marker='.', s=5, zorder=10)
    plt.show()



