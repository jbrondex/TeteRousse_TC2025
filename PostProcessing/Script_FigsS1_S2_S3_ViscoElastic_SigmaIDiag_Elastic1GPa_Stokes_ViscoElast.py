"""
@author: jbrondex

Description:
------------
This file aims at plotting the surface maximum principal stress obtained with the elastic (E=1GPa),
viscous line (A=0.4 MPa-1 a-1), viscous non line (T = Tmap) and visco-elastic frameworks. Each time the background sigma_I
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
    Col_Crevasses_Circ = '#00ff00' ###'#FF0000' ###'#800080' ###Color for representation of circular crevasses
    Col_Crevasses_Other = '#800080' ###Color for representation of non circular crevasses
    Col_Cavity = '#00FFFF'###color for representation of cavity
    Col_Transect = 'lightgrey'###color for representation of transect
    WithCompression = True ###Do we want to plot compressive stresses as well or only tensile stresses ?
    AbsoluteDifference = True ###For timestep sensitivity do we want relative or absolute differente

    if WithCompression:
        Colormap_for_stressmap = 'RdBu_r' ###Colorm map to choose for map of stress
        cmap = cmx.get_cmap(Colormap_for_stressmap )
        clevtight = np.arange(-60.01, 60.01, 5)
        levs_ticktight= np.round(np.arange(-60.01, 60.01, 30))
        clevwide = np.arange(-120.01, 120.01, 0.1)
        levs_tickwide = np.round(np.arange(-120.01,120.01,40))
    else:
        Colormap_for_stressmap = 'viridis' ###Colorm map to choose for map of stress
        cmap = cmx.get_cmap(Colormap_for_stressmap )
        cmap = colors.LinearSegmentedColormap.from_list("viridis_darker", cmap(np.linspace(0.2, 1, 256)))
        clevtight = np.arange(-0.01, 60.01, 5)
        levs_ticktight= np.round(np.arange(0.0,60.01,10))
        clevwide = np.arange(0.0,121,10)
        levs_tickwide = np.round(np.arange(0.0,121,20))

    Colormap_for_comp='BrBG'
    cmapdiff = cmx.get_cmap(Colormap_for_comp )
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
    Col_Names_Elastic = ['Timestep', 'BC', 'NodeNumber', 'X', 'Y', 'Z', 'DispX', 'DispY', 'DispZ', 'vonMises', 'SigmaI','SigmaII','SigmaIII']
    Col_Names_ViscoElastic = ['Timestep', 'Count', 'BC', 'NodeNumber', 'X', 'Y', 'Z', 'DispX', 'DispY', 'DispZ', 'DispVeloX', 'DispVeloY', 'DispVeloZ','vonMises', 'SigmaI','SigmaII','SigmaIII']
    ###We create a single dataframe for each constitutive law : one for elastic, one for viscous line, one for viscous non line
    Data_Simu_El = pd.DataFrame()  ##Initialize an empty dataframe
    Data_Simu_El_P045 = pd.DataFrame()  ##Initialize an empty dataframe
    Data_Simu_ElINCOMP = pd.DataFrame()
    Data_Simu_VL = pd.DataFrame()
    Data_Simu_VNL = pd.DataFrame()
    Data_Simu_LVELCOMP_1s = pd.DataFrame()
    Data_Simu_LVELINCOMP_1s = pd.DataFrame()
    Data_Simu_NLVELCOMP_1s = pd.DataFrame()
    Data_Simu_NLVELINCOMP_1s = pd.DataFrame()
    Data_Simu_VLTRANS_1s = pd.DataFrame()
    Data_Simu_VNLTRANS_1s = pd.DataFrame()
    Data_Simu_LVELCOMP_1mn = pd.DataFrame()
    Data_Simu_LVELINCOMP_1mn = pd.DataFrame()
    Data_Simu_NLVELCOMP_1mn = pd.DataFrame()
    Data_Simu_NLVELINCOMP_1mn = pd.DataFrame()
    Data_Simu_LVELCOMP_5mn = pd.DataFrame()
    Data_Simu_LVELINCOMP_5mn = pd.DataFrame()
    Data_Simu_NLVELCOMP_5mn = pd.DataFrame()
    Data_Simu_NLVELINCOMP_5mn = pd.DataFrame()
    Data_Simu_VLTRANS_5mn = pd.DataFrame()
    Data_Simu_VNLTRANS_5mn = pd.DataFrame()

    for case in ['NoCavity', 'EmptyCavity']:
        filename_el = 'SurfaceOutput_DiagElastic_Poisson03_E1GPa_{}_NoDispLat_v97773a1f2.dat'.format(case)
        filename_el_P045 = 'SurfaceOutput_DiagElastic_Poisson045_E1GPa_{}_NoDispLat_v97773a1f2.dat'.format(case)
        filename_el_incomp = 'SurfaceOutput_DiagElastic_INCOMP_E1GPa_{}_NoDispLat_v97773a1f2.dat'.format(case)

        filename_vl = 'SurfaceOutput_DiagStokes_Line_{}_A04_NoSliding_NoDispLat_v97773a1f2.dat'.format(case)
        filename_vnl = 'SurfaceOutput_DiagStokes_NonLine_{}_Tmap_NoSliding_NoDispLat_v97773a1f2.dat'.format(case)

        filename_lvelcomp_tsp1s = 'SurfaceOutput_DiagLinearViscoElastic_Poisson045_E1GPa_A04_tsp1s_1mn_{}_v97773a1f2_.dat'.format(case)
        filename_lvelincomp_tsp1s = 'SurfaceOutput_DiagLinearViscoElastic_Incomp_E1GPa_A04_tsp1s_5mn_{}_v97773a1f2_.dat'.format(case)
        # # filename_lvelincomp_tsp1s = 'SurfaceOutput_DiagLinearViscoElastic_Incomp_E1GPa_A04_tsp1s_5mn_{}_v97773a1f2_DispMesh_.dat'.format(case)
        filename_nlvelcomp_tsp1s = 'SurfaceOutput_DiagNonLinearViscoElastic_Poisson045_E1GPa_Tmap_tsp1s_1mn_{}_v97773a1f2_.dat'.format(case)
        # filename_nlvelincomp_tsp1s = 'SurfaceOutput_DiagNonLinearViscoElastic_Incomp_E1GPa_Tmap_tsp1s_5mn_{}_v97773a1f2_.dat'.format(case)
        # filename_vl_transient_tsp1s = 'SurfaceOutput_DiagStokesTransient_Line_{}_A04_tsp1s_1mn_v97773a1f2_.dat'.format(case)
        # filename_vnl_transient_tsp1s = 'SurfaceOutput_DiagStokesTransient_NonLine_{}_Tmap_tsp1s_1mn_v97773a1f2_.dat'.format(case)

        # filename_lvelcomp_tsp1mn = 'SurfaceOutput_DiagLinearViscoElastic_Poisson03_E1GPa_A04_tsp1mn_5h_{}_v97773a1f2_.dat'.format(case)
        # filename_lvelincomp_tsp1mn = 'SurfaceOutput_DiagLinearViscoElastic_Incomp_E1GPa_A04_tsp1mn_5h_{}_v97773a1f2_.dat'.format(case)
        # # filename_lvelincomp_tsp1mn = 'SurfaceOutput_DiagLinearViscoElastic_Incomp_E1GPa_A04_tsp1mn_5h_{}_v97773a1f2_DispMesh_.dat'.format(case)
        # filename_nlvelcomp_tsp1mn = 'SurfaceOutput_DiagNonLinearViscoElastic_Poisson03_E1GPa_Tmap_tsp1mn_5h_{}_v97773a1f2_.dat'.format(case)
        # filename_nlvelincomp_tsp1mn = 'SurfaceOutput_DiagNonLinearViscoElastic_Incomp_E1GPa_Tmap_tsp1mn_5h_{}_v97773a1f2_.dat'.format(case)
        #
        filename_lvelcomp_tsp5mn = 'SurfaceOutput_DiagLinearViscoElastic_Poisson045_E1GPa_A04_tsp5mn_1d_{}_v97773a1f2_.dat'.format(case)
        # filename_lvelincomp_tsp5mn = 'SurfaceOutput_DiagLinearViscoElastic_Incomp_E1GPa_A04_tsp5mn_1d_{}_v97773a1f2_.dat'.format(case)
        # # filename_lvelincomp_tsp5mn = 'SurfaceOutput_DiagLinearViscoElastic_Incomp_E1GPa_A04_tsp5mn_1d_{}_v97773a1f2_DispMesh_.dat'.format(case)
        filename_nlvelcomp_tsp5mn = 'SurfaceOutput_DiagNonLinearViscoElastic_Poisson045_E1GPa_Tmap_tsp5mn_1d_{}_v97773a1f2_.dat'.format(case)
        # filename_nlvelincomp_tsp5mn = 'SurfaceOutput_DiagNonLinearViscoElastic_Incomp_E1GPa_Tmap_tsp5mn_1d_{}_v97773a1f2_.dat'.format(case)
        # filename_vl_transient_tsp5mn = 'SurfaceOutput_DiagStokesTransient_Line_{}_A04_tsp5mn_1d_v97773a1f2_.dat'.format(case)
        # filename_vnl_transient_tsp5mn = 'SurfaceOutput_DiagStokesTransient_NonLine_{}_Tmap_tsp5mn_1d_v97773a1f2_.dat'.format(case)

        print('Opening file:', filename_el)
        Data_Simu_El_tmp = pd.read_csv(Pathroot_SimuOutput.joinpath(filename_el), names=Col_Names_Elastic, delim_whitespace=True)
        Data_Simu_El_tmp['Case'] = case
        Data_Simu_El = pd.concat([Data_Simu_El, Data_Simu_El_tmp], ignore_index=True)

        print('Opening file:', filename_el_P045)
        Data_Simu_El_P045_tmp = pd.read_csv(Pathroot_SimuOutput.joinpath(filename_el_P045), names=Col_Names_Elastic, delim_whitespace=True)
        Data_Simu_El_P045_tmp['Case'] = case
        Data_Simu_El_P045 = pd.concat([Data_Simu_El_P045, Data_Simu_El_P045_tmp], ignore_index=True)


        print('Opening file:', filename_el_incomp)
        Data_Simu_ElINCOMP_tmp = pd.read_csv(Pathroot_SimuOutput.joinpath(filename_el_incomp), names=Col_Names_Elastic, delim_whitespace=True)
        Data_Simu_ElINCOMP_tmp['Case'] = case
        Data_Simu_ElINCOMP = pd.concat([Data_Simu_ElINCOMP, Data_Simu_ElINCOMP_tmp], ignore_index=True)

        print('Opening file:', filename_vl)
        Data_Simu_VL_tmp = pd.read_csv(Pathroot_SimuOutput.joinpath(filename_vl), names=Col_Names_Viscous, delim_whitespace=True)
        Data_Simu_VL_tmp['Case'] = case
        Data_Simu_VL = pd.concat([Data_Simu_VL, Data_Simu_VL_tmp], ignore_index=True)

        print('Opening file:', filename_vnl)
        Data_Simu_VNL_tmp = pd.read_csv(Pathroot_SimuOutput.joinpath(filename_vnl), names=Col_Names_Viscous, delim_whitespace=True)
        Data_Simu_VNL_tmp['Case'] = case
        Data_Simu_VNL = pd.concat([Data_Simu_VNL, Data_Simu_VNL_tmp], ignore_index=True)

        print('Opening file:', filename_lvelcomp_tsp1s)
        Data_Simu_LVELCOMP_tmp_1s = pd.read_csv(Pathroot_SimuOutput.joinpath(filename_lvelcomp_tsp1s), names=Col_Names_ViscoElastic, delim_whitespace=True)
        Data_Simu_LVELCOMP_tmp_1s['Case'] = case
        Data_Simu_LVELCOMP_1s = pd.concat([Data_Simu_LVELCOMP_1s, Data_Simu_LVELCOMP_tmp_1s], ignore_index=True)

        print('Opening file:', filename_lvelincomp_tsp1s)
        Data_Simu_LVELINCOMP_tmp_1s = pd.read_csv(Pathroot_SimuOutput.joinpath(filename_lvelincomp_tsp1s), names=Col_Names_ViscoElastic, delim_whitespace=True)
        Data_Simu_LVELINCOMP_tmp_1s['Case'] = case
        Data_Simu_LVELINCOMP_1s = pd.concat([Data_Simu_LVELINCOMP_1s, Data_Simu_LVELINCOMP_tmp_1s], ignore_index=True)

        print('Opening file:', filename_nlvelcomp_tsp1s)
        Data_Simu_NLVELCOMP_tmp_1s = pd.read_csv(Pathroot_SimuOutput.joinpath(filename_nlvelcomp_tsp1s), names=Col_Names_ViscoElastic, delim_whitespace=True)
        Data_Simu_NLVELCOMP_tmp_1s['Case'] = case
        Data_Simu_NLVELCOMP_1s = pd.concat([Data_Simu_NLVELCOMP_1s, Data_Simu_NLVELCOMP_tmp_1s], ignore_index=True)

        # print('Opening file:', filename_nlvelincomp_tsp1s)
        # Data_Simu_NLVELINCOMP_tmp_1s = pd.read_csv(Pathroot_SimuOutput.joinpath(filename_nlvelincomp_tsp1s), names=Col_Names_ViscoElastic, delim_whitespace=True)
        # Data_Simu_NLVELINCOMP_tmp_1s['Case'] = case
        # Data_Simu_NLVELINCOMP_1s = pd.concat([Data_Simu_NLVELINCOMP_1s, Data_Simu_NLVELINCOMP_tmp_1s], ignore_index=True)

        # print('Opening file:', filename_vl_transient_tsp1s)
        # Data_Simu_VLTRANS_tmp_1s = pd.read_csv(Pathroot_SimuOutput.joinpath(filename_vl_transient_tsp1s), names=Col_Names_Viscous, delim_whitespace=True)
        # Data_Simu_VLTRANS_tmp_1s['Case'] = case
        # Data_Simu_VLTRANS_1s = pd.concat([Data_Simu_VLTRANS_1s, Data_Simu_VLTRANS_tmp_1s], ignore_index=True)
        #
        # print('Opening file:', filename_vnl_transient_tsp1s)
        # Data_Simu_VNLTRANS_tmp_1s = pd.read_csv(Pathroot_SimuOutput.joinpath(filename_vnl_transient_tsp1s), names=Col_Names_Viscous, delim_whitespace=True)
        # Data_Simu_VNLTRANS_tmp_1s['Case'] = case
        # Data_Simu_VNLTRANS_1s = pd.concat([Data_Simu_VNLTRANS_1s, Data_Simu_VNLTRANS_tmp_1s], ignore_index=True)

        # print('Opening file:', filename_lvelcomp_tsp1mn)
        # Data_Simu_LVELCOMP_tmp_1mn = pd.read_csv(Pathroot_SimuOutput.joinpath(filename_lvelcomp_tsp1mn), names=Col_Names_ViscoElastic, delim_whitespace=True)
        # Data_Simu_LVELCOMP_tmp_1mn['Case'] = case
        # Data_Simu_LVELCOMP_1mn = pd.concat([Data_Simu_LVELCOMP_1mn, Data_Simu_LVELCOMP_tmp_1mn], ignore_index=True)

        # print('Opening file:', filename_lvelincomp_tsp1mn)
        # Data_Simu_LVELINCOMP_tmp_1mn = pd.read_csv(Pathroot_SimuOutput.joinpath(filename_lvelincomp_tsp1mn), names=Col_Names_ViscoElastic, delim_whitespace=True)
        # Data_Simu_LVELINCOMP_tmp_1mn['Case'] = case
        # Data_Simu_LVELINCOMP_1mn = pd.concat([Data_Simu_LVELINCOMP_1mn, Data_Simu_LVELINCOMP_tmp_1mn],ignore_index=True)

        # print('Opening file:', filename_nlvelcomp_tsp1mn)
        # Data_Simu_NLVELCOMP_tmp_1mn = pd.read_csv(Pathroot_SimuOutput.joinpath(filename_nlvelcomp_tsp1mn), names=Col_Names_ViscoElastic, delim_whitespace=True)
        # Data_Simu_NLVELCOMP_tmp_1mn['Case'] = case
        # Data_Simu_NLVELCOMP_1mn = pd.concat([Data_Simu_NLVELCOMP_1mn, Data_Simu_NLVELCOMP_tmp_1mn], ignore_index=True)
        #
        # print('Opening file:', filename_nlvelincomp_tsp1mn)
        # Data_Simu_NLVELINCOMP_tmp_1mn = pd.read_csv(Pathroot_SimuOutput.joinpath(filename_nlvelincomp_tsp1mn), names=Col_Names_ViscoElastic, delim_whitespace=True)
        # Data_Simu_NLVELINCOMP_tmp_1mn['Case'] = case
        # Data_Simu_NLVELINCOMP_1mn = pd.concat([Data_Simu_NLVELINCOMP_1mn, Data_Simu_NLVELINCOMP_tmp_1mn],ignore_index=True)
        #
        print('Opening file:', filename_lvelcomp_tsp5mn)
        Data_Simu_LVELCOMP_tmp_5mn = pd.read_csv(Pathroot_SimuOutput.joinpath(filename_lvelcomp_tsp5mn), names=Col_Names_ViscoElastic, delim_whitespace=True)
        Data_Simu_LVELCOMP_tmp_5mn['Case'] = case
        Data_Simu_LVELCOMP_5mn = pd.concat([Data_Simu_LVELCOMP_5mn, Data_Simu_LVELCOMP_tmp_5mn], ignore_index=True)

        # print('Opening file:', filename_lvelincomp_tsp5mn)
        # Data_Simu_LVELINCOMP_tmp_5mn = pd.read_csv(Pathroot_SimuOutput.joinpath(filename_lvelincomp_tsp5mn), names=Col_Names_ViscoElastic, delim_whitespace=True)
        # Data_Simu_LVELINCOMP_tmp_5mn['Case'] = case
        # Data_Simu_LVELINCOMP_5mn = pd.concat([Data_Simu_LVELINCOMP_5mn, Data_Simu_LVELINCOMP_tmp_5mn], ignore_index=True)
        #
        print('Opening file:', filename_nlvelcomp_tsp5mn)
        Data_Simu_NLVELCOMP_tmp_5mn = pd.read_csv(Pathroot_SimuOutput.joinpath(filename_nlvelcomp_tsp5mn), names=Col_Names_ViscoElastic, delim_whitespace=True)
        Data_Simu_NLVELCOMP_tmp_5mn['Case'] = case
        Data_Simu_NLVELCOMP_5mn = pd.concat([Data_Simu_NLVELCOMP_5mn, Data_Simu_NLVELCOMP_tmp_5mn], ignore_index=True)


        # print('Opening file:', filename_vl_transient_tsp5mn)
        # Data_Simu_VLTRANS_tmp_5mn = pd.read_csv(Pathroot_SimuOutput.joinpath(filename_vl_transient_tsp5mn), names=Col_Names_Viscous, delim_whitespace=True)
        # Data_Simu_VLTRANS_tmp_5mn['Case'] = case
        # Data_Simu_VLTRANS_5mn = pd.concat([Data_Simu_VLTRANS_5mn, Data_Simu_VLTRANS_tmp_5mn], ignore_index=True)
        #
        # print('Opening file:', filename_vnl_transient_tsp5mn)
        # Data_Simu_VNLTRANS_tmp_5mn = pd.read_csv(Pathroot_SimuOutput.joinpath(filename_vnl_transient_tsp5mn), names=Col_Names_Viscous, delim_whitespace=True)
        # Data_Simu_VNLTRANS_tmp_5mn['Case'] = case
        # Data_Simu_VNLTRANS_5mn = pd.concat([Data_Simu_VNLTRANS_5mn, Data_Simu_VNLTRANS_tmp_5mn], ignore_index=True)

        #
        # print('Opening file:', filename_nlvelincomp_tsp5mn)
        # Data_Simu_NLVELINCOMP_tmp_5mn = pd.read_csv(Pathroot_SimuOutput.joinpath(filename_nlvelincomp_tsp5mn), names=Col_Names_ViscoElastic, delim_whitespace=True)
        # Data_Simu_NLVELINCOMP_tmp_5mn['Case'] = case
        # Data_Simu_NLVELINCOMP_5mn = pd.concat([Data_Simu_NLVELINCOMP_5mn, Data_Simu_NLVELINCOMP_tmp_5mn], ignore_index=True)
    # #####################################################################
    # ###     PLOT MAP OF FIRST MAXI PRINCIPAL STRESS AT THE SURFACE    ###
    # #####################################################################
    ###Create refined regular grid on which all SigmaI for all cases will be interpolated
    x = np.arange(np.floor(np.min(Data_Simu_El['X'])) - 5,np.floor(np.max(Data_Simu_El['X'])) + 5, 0.2)
    y = np.arange(np.floor(np.min(Data_Simu_El['Y'])) - 5, np.floor(np.max(Data_Simu_El['Y'])) + 5, 0.2)
    X, Y = np.meshgrid(x, y)
    # ###Prepare the subplot : fig 1 for linear case, fig 2 for non linear case
    nrows = 3 ### row 1 Stokes line A=0.4, row 2 mu =0.3, row 3 incomp
    ncols = 2  ### col 1 elastic 1GPa mu=0.3, col 2 linear viscous elastic
    fig, axes = plt.subplots(nrows, ncols, figsize=(11, 21), sharex=True, sharey=True, constrained_layout=True)
    # fig2, axes2 = plt.subplots(nrows, ncols, figsize=(11, 21), sharex=True, sharey=True, constrained_layout=True) ###Same as Fig 1 but non linear case
    # fig3, axes3 = plt.subplots(nrows, ncols, figsize=(11, 21), sharex=True, sharey=True, constrained_layout=True) ###Just all version of visco-elastic

    # ####Figure 4 aims at showing evolution of sigmaI depending on tsp for all linear visco_elastic incomp
    # nrows4 = 2 ##sub hour, sub-day
    # ncols4 = 4  ###elastic, tsp 1s, tsp 1mn, tsp 1h // tsp 6h, tsp 12h, tsp 24h, viscous
    # fig4, axes4 = plt.subplots(nrows4, ncols4, figsize=(21, 10), sharex=True, sharey=True, constrained_layout=True)

    # ####Figure 4 bis :idem Figure 4 for all non_linear visco_elastic incomp
    # nrows4bis = 2 ##sub hour, sub-day
    # ncols4bis = 4  ###elastic, tsp 1s, tsp 1mn, tsp 1h // tsp 6h, tsp 12h, tsp 24h, viscous
    # fig4bis, axes4bis = plt.subplots(nrows4bis, ncols4bis, figsize=(21, 10), sharex=True, sharey=True, constrained_layout=True)

    ####Figure 5 aims at showing evolution of sigmaI depending on tsp for all linear visco_elastic compressible
    nrows5 = 2 ##sub hour, sub-day
    ncols5 = 4  ###elastic, tsp 1s, tsp 1mn, tsp 1h // tsp 6h, tsp 12h, tsp 24h, none
    fig5, axes5 = plt.subplots(nrows5, ncols5, figsize=(21, 10), sharex=True, sharey=True, constrained_layout=True)

    ####Figure 5 bis : idem Figure 5 for all non_linear visco_elastic compressible
    nrows5bis = 2 ##sub hour, sub-day
    ncols5bis = 4  ###elastic, tsp 1s, tsp 1mn, tsp 1h // tsp 6h, tsp 12h, tsp 24h, none
    fig5bis, axes5bis = plt.subplots(nrows5bis, ncols5bis, figsize=(18, 13), sharex=True, sharey=True, constrained_layout=True)

    # ####Figure 6 aims at showing sensitivity of visco-elastic to timestep : result @1mn with tsp=1s vs tsp=1mn
    # nrows6 = 4 ##Empty cav only but Line comp/incomp and Non Line comp/incomp
    # ncols6 = 3  ###tsp 1s, tsp 1mn, difference
    # fig6, axes6 = plt.subplots(nrows6, ncols6, figsize=(13.5, 24.5), sharex=True, sharey=True, constrained_layout=True)
    #
    # ####Figure 6 bis aims at showing sensitivity of visco-elastic to timestep : result @5mn with tsp=1mn vs tsp=5mn
    # nrows6bis = 4 ##Empty cav only but Line comp/incomp and Non Line comp/incomp
    # ncols6bis = 3  ###tsp 1s, tsp 1mn, difference
    # fig6bis, axes6bis = plt.subplots(nrows6bis, ncols6bis, figsize=(13.5, 24.5), sharex=True, sharey=True, constrained_layout=True)
    #
    # ####Figure 7 aims at showing sensitivity of visco-elastic to timestep : result @1h with tsp=1mn vs tsp=10mn
    # nrows7= 4 ##Empty cav only but Line comp/incomp and Non Line comp/incomp
    # ncols7 = 3  ###tsp 1s, tsp 1mn, difference
    # fig7, axes7 = plt.subplots(nrows7, ncols7, figsize=(13.5, 24.5), sharex=True, sharey=True, constrained_layout=True)

    # ####Figure 8 aims at showing evolution of sigmaI depending on tsp for all linear viscous transient
    # nrows8 = 2 ##sub hour, sub-day
    # ncols8 = 4  ###elastic, tsp 1s, tsp 1mn, tsp 1h // tsp 6h, tsp 12h, tsp 24h, none
    # fig8, axes8 = plt.subplots(nrows8, ncols8, figsize=(21, 10), sharex=True, sharey=True, constrained_layout=True)
    #
    # ####Figure 8 bis : idem Figure 8 for all non_linear glen transient
    # nrows8bis = 2 ##sub hour, sub-day
    # ncols8bis = 4  ###elastic, tsp 1s, tsp 1mn, tsp 1h // tsp 6h, tsp 12h, tsp 24h, none
    # fig8bis, axes8bis = plt.subplots(nrows8bis, ncols8bis, figsize=(18, 13), sharex=True, sharey=True, constrained_layout=True)

    ###### START PLOT FIGURE 1 ######
    fig.suptitle('Linear Viscosity', fontsize=23, weight='bold')
    for k,law in enumerate(['ElP045','LVELCOMP', 'ElINCOMP', 'LVELINCOMP','none','VL']):
        df_cav =pd.DataFrame()
        df_nocav=pd.DataFrame()
        if law == 'ElP045':
            df_nocav =  Data_Simu_El_P045[Data_Simu_El_P045['Case']=='NoCavity']
            df_cav = Data_Simu_El_P045[Data_Simu_El_P045['Case']=='EmptyCavity']
            Title = r'Elastic $\nu = 0.45$'
        elif law == 'ElINCOMP':
            df_nocav =  Data_Simu_ElINCOMP[Data_Simu_ElINCOMP['Case']=='NoCavity']
            df_cav = Data_Simu_ElINCOMP[Data_Simu_ElINCOMP['Case']=='EmptyCavity']
            Title = r'Elastic incomp.'
        elif law == 'VL':
            df_nocav =  Data_Simu_VL[Data_Simu_VL['Case']=='NoCavity']
            df_cav = Data_Simu_VL[Data_Simu_VL['Case']=='EmptyCavity']
            Title = 'Pure Viscous incomp.'
        elif law=='LVELCOMP':
            df_nocav =  Data_Simu_LVELCOMP_1s[(Data_Simu_LVELCOMP_1s['Timestep']==1) & (Data_Simu_LVELCOMP_1s['Case']=='NoCavity')]
            df_cav = Data_Simu_LVELCOMP_1s[(Data_Simu_LVELCOMP_1s['Timestep']==1) & (Data_Simu_LVELCOMP_1s['Case']=='EmptyCavity')]
            Title = r'Lin. Viscoelastic $\nu = 0.45$ @1s'
        elif law=='LVELINCOMP':
            df_nocav =  Data_Simu_LVELINCOMP_1s[(Data_Simu_LVELINCOMP_1s['Timestep']==1) & (Data_Simu_LVELINCOMP_1s['Case']=='NoCavity')]
            df_cav = Data_Simu_LVELINCOMP_1s[(Data_Simu_LVELINCOMP_1s['Timestep']==1) & (Data_Simu_LVELINCOMP_1s['Case']=='EmptyCavity')]
            Title = r'Lin. Viscoelastic incomp. @1s'
        ##get proper ax
        i, j = divmod(k, 2)
        ax = axes[i, j]
        if k == 4:
            ax.axis("off")
            continue
        if ((i == 2 and j==1) or (i == 1 and j==0)):
            ax.set_xlabel(r'X [km]', fontsize=22)
        if ((i==2 and j==1) or j == 0):
            ax.set_ylabel(r'Y [km]', fontsize=22)
        if i==0 and j == 0:
            ax.text(-0.45, 0.5, r'Compressible', transform=ax.transAxes, rotation=90, ha='center', va='center', fontsize=18, weight='bold')
        if i==1 and j == 0:
            ax.text(-0.45, 0.5, r'Incompressible', transform=ax.transAxes, rotation=90, ha='center', va='center', fontsize=18, weight='bold')
        ax.tick_params(labelsize=18)  # fontsize of the tick labels
        ax.grid(True)
        ax.grid(alpha=0.5)
        ax.set_title(Title, fontsize=16)
        ###Force x and y axis to same scale
        ax.set_aspect('equal', adjustable='box')
        ###SigmaI is max principal stress and at the surface always one of the principal stress is zero
        ### Therefore SigmaI is either positive or zero. Negative values of SigmaI are numerical artefacts
        df_nocav.loc[df_nocav["SigmaI"] < 0, "SigmaI"] = 0
        df_cav.loc[df_cav["SigmaI"] < 0, "SigmaI"] = 0
        ###Poject SigmaI on interpolation grid
        SigmaI_NoCav = Interpolate_field(df_nocav, 'SigmaI', X, Y)
        SigmaI_Cav = Interpolate_field(df_cav, 'SigmaI', X, Y)
        ###shading
        clevs=clevtight ## cbar for shading
        #colorbar
        levs_ticks=levs_ticktight
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
                continue ##We don't plot non circular crevasses anymre. Before we had: col = Col_Crevasses_Other
            ax.plot(Df_plot['X'].values/1000, Df_plot['Y'].values/1000, color=col, linestyle='-', linewidth=2)
        ###Plot cavity contour
        ax.plot(cavity_contour[:, 0]/1000,cavity_contour[:, 1]/1000,color=Col_Cavity,linestyle='-',linewidth=2.5)
    ###Make a horizontal color bar common to all subplots
    cbar = fig.colorbar(CS1,
                        ax=axes.ravel().tolist(),  # all subplots
                        orientation='horizontal',
                        ticks=levs_ticks,
                        fraction=0.04,
                        pad=0.08)
    cbar.set_label(r'$\Delta \sigma_\mathrm{I}$ [kPa]', fontsize=20)
    cbar.ax.tick_params(labelsize=18)  # Adjust tick label size

    # ###### START PLOT FIGURE 2 ######
    # fig2.suptitle('Glen\'s Viscosity', fontsize=23, weight='bold')
    # for k,law in enumerate(['none','VNL','El','NLVELCOMP', 'ElINCOMP', 'NLVELINCOMP']):
    #     df_cav2 =pd.DataFrame()
    #     df_nocav2=pd.DataFrame()
    #     if law == 'El':
    #         df_nocav2 =  Data_Simu_El[Data_Simu_El['Case']=='NoCavity']
    #         df_cav2 = Data_Simu_El[Data_Simu_El['Case']=='EmptyCavity']
    #         Title = r'Elastic $\nu = 0.3$'
    #     elif law == 'ElINCOMP':
    #         df_nocav2 =  Data_Simu_ElINCOMP[Data_Simu_ElINCOMP['Case']=='NoCavity']
    #         df_cav2 = Data_Simu_ElINCOMP[Data_Simu_ElINCOMP['Case']=='EmptyCavity']
    #         Title = r'Elastic incomp.'
    #     elif law == 'VNL':
    #         df_nocav2 =  Data_Simu_VNL[Data_Simu_VNL['Case']=='NoCavity']
    #         df_cav2 = Data_Simu_VNL[Data_Simu_VNL['Case']=='EmptyCavity']
    #         Title = 'Glen'
    #     elif law=='NLVELCOMP':
    #         df_nocav2 =  Data_Simu_NLVELCOMP_1s[(Data_Simu_NLVELCOMP_1s['Timestep']==1) & (Data_Simu_NLVELCOMP_1s['Case']=='NoCavity')]
    #         df_cav2 = Data_Simu_NLVELCOMP_1s[(Data_Simu_NLVELCOMP_1s['Timestep']==1) & (Data_Simu_NLVELCOMP_1s['Case']=='EmptyCavity')]
    #         Title = r'NL Visco-el $\nu = 0.3$ @1s'
    #     elif law=='NLVELINCOMP':
    #         df_nocav2 =  Data_Simu_NLVELINCOMP_1s[(Data_Simu_NLVELINCOMP_1s['Timestep']==1) & (Data_Simu_NLVELINCOMP_1s['Case']=='NoCavity')]
    #         df_cav2 = Data_Simu_NLVELINCOMP_1s[(Data_Simu_NLVELINCOMP_1s['Timestep']==1) & (Data_Simu_NLVELINCOMP_1s['Case']=='EmptyCavity')]
    #         Title = r'NL Visco-el incomp. @1s'
    #     ##get proper ax
    #     i, j = divmod(k, 2)
    #     ax2 = axes2[i, j]
    #     if k == 0:
    #         ax2.axis("off")
    #         continue
    #     if i == 2:
    #         ax2.set_xlabel(r'X [km]', fontsize=22)
    #     if ((i == 0 and j == 1) or j == 0):
    #         ax2.set_ylabel(r'Y [km]', fontsize=22)
    #     if i == 1 and j == 0:
    #         ax2.text(-0.45, 0.5, r'Compressible', transform=ax2.transAxes, rotation=90, ha='center', va='center',
    #                 fontsize=18, weight='bold')
    #     if i == 2 and j == 0:
    #         ax2.text(-0.45, 0.5, r'Incompressible', transform=ax2.transAxes, rotation=90, ha='center', va='center',
    #                 fontsize=18, weight='bold')
    #     ax2.tick_params(labelsize=18)  # fontsize of the tick labels
    #     ax2.grid(True)
    #     ax2.grid(alpha=0.5)
    #     ax2.set_title(Title, fontsize=16)
    #     ###Force x and y axis to same scale
    #     ax2.set_aspect('equal', adjustable='box')
    #     ###Poject SigmaI on interpolation grid
    #     SigmaI_NoCav = Interpolate_field(df_nocav2, 'SigmaI', X, Y)
    #     SigmaI_Cav = Interpolate_field(df_cav2, 'SigmaI', X, Y)
    #     ###shading
    #     clevs2=np.arange(-0.01,81,5) ## cbar for shading
    #     #colorbar
    #     levs_ticks2=np.arange(0.0,81,10)
    #     ###Fills up the map with colors for SigmaEq
    #     CS2 = ax2.contourf(X/1000, Y/1000, (SigmaI_Cav-SigmaI_NoCav)*1000, clevs2, cmap=cmap,extend='both')
    #     ax2.plot(xc / 1000, yc / 1000, color='k', linewidth=2)
    #     ###Below we remove colors that are outside of glacier contour
    #     clippath = mpltPath(np.c_[xc/1000, yc/1000])
    #     patch = PathPatch(clippath, facecolor='none')
    #     # ax = plt.gca()
    #     ax2.add_patch(patch)
    #     for c in [CS2]:
    #         c.set_clip_path(patch)
    #     ##plot crevasses as continuous line
    #     for crev_num in np.arange(Df_Crevasses['Crevasse Number'].min(), Df_Crevasses['Crevasse Number'].max() + 1):
    #         ###get proper points
    #         Df_plot = Df_Crevasses[Df_Crevasses['Crevasse Number'] == crev_num]
    #         if Df_plot['IsCircular'].all():  ##different colors for circular crevasses and other crevasses
    #             col = Col_Crevasses_Circ
    #         else:
    #             continue ##We don't plot non circular crevasses anymre. Before we had: col = Col_Crevasses_Other
    #         ax2.plot(Df_plot['X'].values/1000, Df_plot['Y'].values/1000, color=col, linestyle='-', linewidth=2)
    #     ###Plot cavity contour
    #     ax2.plot(cavity_contour[:, 0]/1000,cavity_contour[:, 1]/1000,color=Col_Cavity,linestyle='-',linewidth=2.5)
    # ###Make a horizontal color bar common to all subplots
    # cbar2 = fig2.colorbar(CS2,
    #                     ax=axes2.ravel().tolist(),  # all subplots
    #                     orientation='horizontal',
    #                     ticks=levs_ticks2,
    #                     fraction=0.04,
    #                     pad=0.08)
    # cbar2.set_label(r'$\sigma_\mathrm{I}$ [kPa]', fontsize=20)
    # cbar2.ax.tick_params(labelsize=18)  # Adjust tick label size
    #
    # ###### START PLOT FIGURE 3 ######
    # fig3.suptitle('Linear vs. Non Linear Visco-elastic @1s', fontsize=23, weight='bold')
    # for k,law in enumerate(['LVELCOMP', 'LVELINCOMP','NLVELCOMP', 'NLVELINCOMP', 'El', 'ElINCOMP']):
    #     df_cav3 =pd.DataFrame()
    #     df_nocav3=pd.DataFrame()
    #     if law == 'LVELCOMP':
    #         df_nocav3 =  Data_Simu_LVELCOMP_1s[(Data_Simu_LVELCOMP_1s['Timestep']==1) & (Data_Simu_LVELCOMP_1s['Case']=='NoCavity')]
    #         df_cav3 = Data_Simu_LVELCOMP_1s[(Data_Simu_LVELCOMP_1s['Timestep']==1) & (Data_Simu_LVELCOMP_1s['Case']=='EmptyCavity')]
    #         Title = r'Linear Visco-elastic $\nu = 0.3$'
    #     elif law == 'LVELINCOMP':
    #         df_nocav3 =  Data_Simu_LVELINCOMP_1s[(Data_Simu_LVELINCOMP_1s['Timestep']==1) & (Data_Simu_LVELINCOMP_1s['Case']=='NoCavity')]
    #         df_cav3 = Data_Simu_LVELINCOMP_1s[(Data_Simu_LVELINCOMP_1s['Timestep']==1) & (Data_Simu_LVELINCOMP_1s['Case']=='EmptyCavity')]
    #         Title = r'Linear Visco-elastic incomp.'
    #     elif law=='NLVELCOMP':
    #         df_nocav3 =  Data_Simu_NLVELCOMP_1s[(Data_Simu_NLVELCOMP_1s['Timestep']==1) & (Data_Simu_NLVELCOMP_1s['Case']=='NoCavity')]
    #         df_cav3 = Data_Simu_NLVELCOMP_1s[(Data_Simu_NLVELCOMP_1s['Timestep']==1) & (Data_Simu_NLVELCOMP_1s['Case']=='EmptyCavity')]
    #         Title = r'NL Visco-elastic $\nu = 0.3$'
    #     elif law=='NLVELINCOMP':
    #         df_nocav3 =  Data_Simu_NLVELINCOMP_1s[(Data_Simu_NLVELINCOMP_1s['Timestep']==1) & (Data_Simu_NLVELINCOMP_1s['Case']=='NoCavity')]
    #         df_cav3 = Data_Simu_NLVELINCOMP_1s[(Data_Simu_NLVELINCOMP_1s['Timestep']==1) & (Data_Simu_NLVELINCOMP_1s['Case']=='EmptyCavity')]
    #         Title = r'NL Visco-elastic incomp.'
    #     elif law == 'El':
    #         df_nocav3 =  Data_Simu_El[Data_Simu_El['Case']=='NoCavity']
    #         df_cav3 = Data_Simu_El[Data_Simu_El['Case']=='EmptyCavity']
    #         Title = r'Elastic $\nu = 0.3$'
    #     elif law == 'ElINCOMP':
    #         df_nocav3 =  Data_Simu_ElINCOMP[Data_Simu_ElINCOMP['Case']=='NoCavity']
    #         df_cav3 = Data_Simu_ElINCOMP[Data_Simu_ElINCOMP['Case']=='EmptyCavity']
    #         Title = r'Elastic incomp.'
    #     ##get proper ax
    #     i, j = divmod(k, 2)
    #     ax3 = axes3[i, j]
    #     if i == 2:
    #         ax3.set_xlabel(r'X [km]', fontsize=22)
    #     if j == 0:
    #         ax3.set_ylabel(r'Y [km]', fontsize=22)
    #     ax3.tick_params(labelsize=18)  # fontsize of the tick labels
    #     ax3.grid(alpha=0.5)
    #     ax3.set_title(Title, fontsize=16)
    #     ###Force x and y axis to same scale
    #     ax3.set_aspect('equal', adjustable='box')
    #     ###Poject SigmaI on interpolation grid
    #     SigmaI_NoCav = Interpolate_field(df_nocav3, 'SigmaI', X, Y)
    #     SigmaI_Cav = Interpolate_field(df_cav3, 'SigmaI', X, Y)
    #     ###shading
    #     clevs3=np.arange(-0.01,51,5) ## cbar for shading
    #     #colorbar
    #     levs_ticks3=np.arange(0.0,51,10)
    #     ###Fills up the map with colors for SigmaEq
    #     CS3 = ax3.contourf(X/1000, Y/1000, (SigmaI_Cav-SigmaI_NoCav)*1000, clevs3, cmap=cmap,extend='both')
    #     ax3.plot(xc / 1000, yc / 1000, color='k', linewidth=2)
    #     ###Below we remove colors that are outside of glacier contour
    #     clippath = mpltPath(np.c_[xc/1000, yc/1000])
    #     patch = PathPatch(clippath, facecolor='none')
    #     # ax = plt.gca()
    #     ax3.add_patch(patch)
    #     for c in [CS3]:
    #         c.set_clip_path(patch)
    #     ##plot crevasses as continuous line
    #     for crev_num in np.arange(Df_Crevasses['Crevasse Number'].min(), Df_Crevasses['Crevasse Number'].max() + 1):
    #         ###get proper points
    #         Df_plot = Df_Crevasses[Df_Crevasses['Crevasse Number'] == crev_num]
    #         if Df_plot['IsCircular'].all():  ##different colors for circular crevasses and other crevasses
    #             col = Col_Crevasses_Circ
    #         else:
    #             continue ##We don't plot non circular crevasses anymre. Before we had: col = Col_Crevasses_Other
    #         ax3.plot(Df_plot['X'].values/1000, Df_plot['Y'].values/1000, color=col, linestyle='-', linewidth=2)
    #     ###Plot cavity contour
    #     ax3.plot(cavity_contour[:, 0]/1000,cavity_contour[:, 1]/1000,color=Col_Cavity,linestyle='-',linewidth=2.5)
    #     ###Make a horizontal color bar common to all subplots
    # cbar3 = fig3.colorbar(CS3,
    #                     ax=axes3.ravel().tolist(),  # all subplots
    #                     orientation='horizontal',
    #                     ticks=levs_ticks3,
    #                     fraction=0.04,
    #                     pad=0.08)
    # cbar3.set_label(r'$\sigma_\mathrm{I}$ [kPa]', fontsize=20)
    # cbar3.ax.tick_params(labelsize=18)  # Adjust tick label size
    #

    # ###### START PLOT FIGURE 4 ######
    # fig4.suptitle('Incompressible Linear Visco-elastic with Mesh Displacement', fontsize=23, weight='bold')
    # print('starting Figure 4')
    # for i,law in enumerate(['Subhour', 'Subday']):
    #     for j,tsp in enumerate(['a','b','c','d']):
    #         ##get proper ax
    #         ax4 = axes4[i,j]
    #         ###Get data for each subplots one by one
    #         df_cav4 =pd.DataFrame()
    #         df_nocav4=pd.DataFrame()
    #         if i==0:
    #             if j == 0:##Linear elastic
    #                 ax4.set_ylabel(r'Y [km]', fontsize=22)
    #                 df_nocav4 = Data_Simu_ElINCOMP[Data_Simu_ElINCOMP['Case'] == 'NoCavity']
    #                 df_cav4 = Data_Simu_ElINCOMP[Data_Simu_ElINCOMP['Case'] == 'EmptyCavity']
    #                 Title = r'Elastic incomp.'
    #             elif j == 1:##Linear visco elastic after 1s
    #                 df_nocav4 =  Data_Simu_LVELINCOMP_1s[(Data_Simu_LVELINCOMP_1s['Timestep']==1) & (Data_Simu_LVELINCOMP_1s['Case']=='NoCavity')]
    #                 df_cav4 = Data_Simu_LVELINCOMP_1s[(Data_Simu_LVELINCOMP_1s['Timestep']==1) & (Data_Simu_LVELINCOMP_1s['Case']=='EmptyCavity')]
    #                 Title = r'@1s tsp = 1s'
    #             elif j == 2: ##Linear visco elastic after 1mn
    #                 df_nocav4 =  Data_Simu_LVELINCOMP_1s[(Data_Simu_LVELINCOMP_1s['Timestep']==60) & (Data_Simu_LVELINCOMP_1s['Case']=='NoCavity')]
    #                 df_cav4 = Data_Simu_LVELINCOMP_1s[(Data_Simu_LVELINCOMP_1s['Timestep']==60) & (Data_Simu_LVELINCOMP_1s['Case']=='EmptyCavity')]
    #                 Title = r'@1mn tsp = 1s'
    #             elif j == 3: ##Linear visco elastic after 1h
    #                 df_nocav4 =  Data_Simu_LVELINCOMP_1mn[(Data_Simu_LVELINCOMP_1mn['Timestep']==60) & (Data_Simu_LVELINCOMP_1mn['Case']=='NoCavity')]
    #                 df_cav4 = Data_Simu_LVELINCOMP_1mn[(Data_Simu_LVELINCOMP_1mn['Timestep']==60) & (Data_Simu_LVELINCOMP_1mn['Case']=='EmptyCavity')]
    #                 Title = r'@1h tsp = 1mn'
    #         elif i==1:
    #             ax4.set_xlabel(r'X [km]', fontsize=22)
    #             if j == 0:##Linear visco elastic after 6h
    #                 ax4.set_ylabel(r'Y [km]', fontsize=22)
    #                 df_nocav4 =  Data_Simu_LVELINCOMP_5mn[(Data_Simu_LVELINCOMP_5mn['Timestep']==6*60/5) & (Data_Simu_LVELINCOMP_5mn['Case']=='NoCavity')]
    #                 df_cav4 = Data_Simu_LVELINCOMP_5mn[(Data_Simu_LVELINCOMP_5mn['Timestep']==6*60/5) & (Data_Simu_LVELINCOMP_5mn['Case']=='EmptyCavity')]
    #                 Title = r'@6h tsp = 5mn'
    #             elif j == 1:##Linear visco elastic tsp = 12h
    #                 df_nocav4 =  Data_Simu_LVELINCOMP_5mn[(Data_Simu_LVELINCOMP_5mn['Timestep']==12*60/5) & (Data_Simu_LVELINCOMP_5mn['Case']=='NoCavity')]
    #                 df_cav4 = Data_Simu_LVELINCOMP_5mn[(Data_Simu_LVELINCOMP_5mn['Timestep']==12*60/5) & (Data_Simu_LVELINCOMP_5mn['Case']=='EmptyCavity')]
    #                 Title = r'@12h tsp = 5mn'
    #             elif j == 2: ##Linear visco elastic tsp = 24h
    #                 df_nocav4 =  Data_Simu_LVELINCOMP_5mn[(Data_Simu_LVELINCOMP_5mn['Timestep']==24*60/5) & (Data_Simu_LVELINCOMP_5mn['Case']=='NoCavity')]
    #                 df_cav4 = Data_Simu_LVELINCOMP_5mn[(Data_Simu_LVELINCOMP_5mn['Timestep']==24*60/5) & (Data_Simu_LVELINCOMP_5mn['Case']=='EmptyCavity')]
    #                 Title = r'@24h tsp = 5mn'
    #             elif j == 3: ##pure linear viscous incomp
    #                 df_nocav4 = Data_Simu_VL[Data_Simu_VL['Case'] == 'NoCavity']
    #                 df_cav4 = Data_Simu_VL[Data_Simu_VL['Case'] == 'EmptyCavity']
    #                 Title = 'L. Viscous incomp.'
    #         ###Fill up the plots
    #         ax4.set_title(Title, fontsize=18, weight='bold')
    #         ax4.tick_params(labelsize=18)  # fontsize of the tick labels
    #         ax4.grid(True)
    #         ax4.grid(alpha=0.5)
    #         ###Force x and y axis to same scale
    #         ax4.set_aspect('equal', adjustable='box')
    #         ###Poject SigmaI on interpolation grid
    #         SigmaI_NoCav = Interpolate_field(df_nocav4, 'SigmaI', X, Y)
    #         SigmaI_Cav = Interpolate_field(df_cav4, 'SigmaI', X, Y)
    #         ###shading
    #         clevs4=clevtight ## cbar for shading
    #         #colorbar
    #         levs_ticks4=levs_ticktight
    #         ###Fills up the map with colors for SigmaEq
    #         CS4= ax4.contourf(X/1000, Y/1000, (SigmaI_Cav-SigmaI_NoCav)*1000, clevs4, cmap=cmap,extend='both')
    #         ax4.plot(xc / 1000, yc / 1000, color='k', linewidth=2)
    #
    #         ###Below we remove colors that are outside of glacier contour
    #         clippath = mpltPath(np.c_[xc/1000, yc/1000])
    #         patch = PathPatch(clippath, facecolor='none')
    #         # ax = plt.gca()
    #         ax4.add_patch(patch)
    #         for c in [CS4]:
    #             c.set_clip_path(patch)
    #         ##plot crevasses as continuous line
    #         for crev_num in np.arange(Df_Crevasses['Crevasse Number'].min(), Df_Crevasses['Crevasse Number'].max() + 1):
    #             ###get proper points
    #             Df_plot = Df_Crevasses[Df_Crevasses['Crevasse Number'] == crev_num]
    #             if Df_plot['IsCircular'].all():  ##different colors for circular crevasses and other crevasses
    #                 col = Col_Crevasses_Circ
    #             else:
    #                 continue ##We don't plot non circular crevasses anymre. Before we had: col = Col_Crevasses_Other
    #             ax4.plot(Df_plot['X'].values/1000, Df_plot['Y'].values/1000, color=col, linestyle='-', linewidth=2)
    #         ###Plot cavity contour
    #         ax4.plot(cavity_contour[:, 0]/1000,cavity_contour[:, 1]/1000,color=Col_Cavity,linestyle='-',linewidth=2.5)
    #         ###Make a horizontal color bar common to all subplots
    # cbar4 = fig4.colorbar(CS4,
    #                     ax=axes4.ravel().tolist(),  # all subplots
    #                     orientation='horizontal',
    #                     ticks=levs_ticks4,
    #                     fraction=0.04,
    #                     pad=0.08)
    # cbar4.set_label(r'$\sigma_\mathrm{I}$ [kPa]', fontsize=20)
    # cbar4.ax.tick_params(labelsize=18)  # Adjust tick label size

    # ###### START PLOT FIGURE 4 BIS ######
    # fig4bis.suptitle('Incompressible Non-Linear Visco-elastic', fontsize=23, weight='bold')
    # print('starting Figure 4 bis')
    # for i,law in enumerate(['Subhour', 'Subday']):
    #     for j,tsp in enumerate(['a','b','c','d']):
    #         ##get proper ax
    #         ax4 = axes4bis[i,j]
    #         ###Get data for each subplots one by one
    #         df_cav4 =pd.DataFrame()
    #         df_nocav4=pd.DataFrame()
    #         if i==0:
    #             if j == 0:##Linear elastic
    #                 ax4.set_ylabel(r'Y [km]', fontsize=22)
    #                 df_nocav4 = Data_Simu_ElINCOMP[Data_Simu_ElINCOMP['Case'] == 'NoCavity']
    #                 df_cav4 = Data_Simu_ElINCOMP[Data_Simu_ElINCOMP['Case'] == 'EmptyCavity']
    #                 Title = r'Elastic incomp.'
    #             elif j == 1:##Linear visco elastic after 1s tsp =1s
    #                 df_nocav4 =  Data_Simu_NLVELINCOMP_1s[(Data_Simu_NLVELINCOMP_1s['Timestep']==1) & (Data_Simu_NLVELINCOMP_1s['Case']=='NoCavity')]
    #                 df_cav4 = Data_Simu_NLVELINCOMP_1s[(Data_Simu_NLVELINCOMP_1s['Timestep']==1) & (Data_Simu_NLVELINCOMP_1s['Case']=='EmptyCavity')]
    #                 Title = r'@1s tsp = 1s'
    #             elif j == 2: ##Linear visco elastic after 1mn tsp = 1s
    #                 df_nocav4 =  Data_Simu_NLVELINCOMP_1s[(Data_Simu_NLVELINCOMP_1s['Timestep']==60) & (Data_Simu_NLVELINCOMP_1s['Case']=='NoCavity')]
    #                 df_cav4 = Data_Simu_NLVELINCOMP_1s[(Data_Simu_NLVELINCOMP_1s['Timestep']==60) & (Data_Simu_NLVELINCOMP_1s['Case']=='EmptyCavity')]
    #                 Title = r'@1mn tsp = 1s'
    #             elif j == 3: ##Linear visco elastic after 1h tsp = 1mn
    #                 df_nocav4 =  Data_Simu_NLVELINCOMP_1mn[(Data_Simu_NLVELINCOMP_1mn['Timestep']==60) & (Data_Simu_NLVELINCOMP_1mn['Case']=='NoCavity')]
    #                 df_cav4 = Data_Simu_NLVELINCOMP_1mn[(Data_Simu_NLVELINCOMP_1mn['Timestep']==60) & (Data_Simu_NLVELINCOMP_1mn['Case']=='EmptyCavity')]
    #                 Title = r'@1h tsp = 1mn'
    #         elif i==1:
    #             ax4.set_xlabel(r'X [km]', fontsize=22)
    #             if j == 0:##Linear visco elastic after 6h tsp =5mn
    #                 ax4.set_ylabel(r'Y [km]', fontsize=22)
    #                 df_nocav4 = Data_Simu_NLVELINCOMP_5mn[(Data_Simu_NLVELINCOMP_5mn['Timestep'] == 6 * 60 / 5) & (Data_Simu_NLVELINCOMP_5mn['Case'] == 'NoCavity')]
    #                 df_cav4 = Data_Simu_NLVELINCOMP_5mn[(Data_Simu_NLVELINCOMP_5mn['Timestep'] == 6 * 60 / 5) & (Data_Simu_NLVELINCOMP_5mn['Case'] == 'EmptyCavity')]
    #                 Title = r'@6h tsp = 5mn'
    #             elif j == 1:  ##Linear visco elastic tsp = 12h
    #                 df_nocav4 = Data_Simu_NLVELINCOMP_5mn[(Data_Simu_NLVELINCOMP_5mn['Timestep'] == 12 * 60 / 5) & (Data_Simu_NLVELINCOMP_5mn['Case'] == 'NoCavity')]
    #                 df_cav4 = Data_Simu_NLVELINCOMP_5mn[(Data_Simu_NLVELINCOMP_5mn['Timestep'] == 12 * 60 / 5) & (Data_Simu_NLVELINCOMP_5mn['Case'] == 'EmptyCavity')]
    #                 Title = r'@12h tsp = 5mn'
    #             elif j == 2:  ##Linear visco elastic tsp = 24h
    #                 df_nocav4 = Data_Simu_NLVELINCOMP_5mn[(Data_Simu_NLVELINCOMP_5mn['Timestep'] == 24 * 60 / 5) & (Data_Simu_NLVELINCOMP_5mn['Case'] == 'NoCavity')]
    #                 df_cav4 = Data_Simu_NLVELINCOMP_5mn[(Data_Simu_NLVELINCOMP_5mn['Timestep'] == 24 * 60 / 5) & (Data_Simu_NLVELINCOMP_5mn['Case'] == 'EmptyCavity')]
    #                 Title = r'@24h tsp = 5mn'
    #             elif j == 3: ##pure linear viscous incomp
    #                 df_nocav4 = Data_Simu_VNL[Data_Simu_VNL['Case'] == 'NoCavity']
    #                 df_cav4 = Data_Simu_VNL[Data_Simu_VNL['Case'] == 'EmptyCavity']
    #                 Title = 'Glen'
    #         ###Fill up the plots
    #         ax4.set_title(Title, fontsize=18, weight='bold')
    #         ax4.tick_params(labelsize=18)  # fontsize of the tick labels
    #         ax4.grid(True)
    #         ax4.grid(alpha=0.5)
    #         ###Force x and y axis to same scale
    #         ax4.set_aspect('equal', adjustable='box')
    #         ###Poject SigmaI on interpolation grid
    #         ###Fills up the map with colors for SigmaEq
    #         SigmaI_NoCav = Interpolate_field(df_nocav4, 'SigmaI', X, Y)
    #         SigmaI_Cav = Interpolate_field(df_cav4, 'SigmaI', X, Y)
    #         ###shading
    #         clevs4=clevwide ## cbar for shading
    #         #colorbar
    #         levs_ticks4=levs_tickwide
    #         ###Fills up the map with colors for SigmaEq
    #         CS4= ax4.contourf(X/1000, Y/1000, (SigmaI_Cav-SigmaI_NoCav)*1000, clevs4, cmap=cmap,extend='both')
    #         ax4.plot(xc / 1000, yc / 1000, color='k', linewidth=2)
    #         ###Below we remove colors that are outside of glacier contour
    #         clippath = mpltPath(np.c_[xc/1000, yc/1000])
    #         patch = PathPatch(clippath, facecolor='none')
    #         # ax = plt.gca()
    #         ax4.add_patch(patch)
    #         for c in [CS4]:
    #             c.set_clip_path(patch)
    #         ##plot crevasses as continuous line
    #         for crev_num in np.arange(Df_Crevasses['Crevasse Number'].min(), Df_Crevasses['Crevasse Number'].max() + 1):
    #             ###get proper points
    #             Df_plot = Df_Crevasses[Df_Crevasses['Crevasse Number'] == crev_num]
    #             if Df_plot['IsCircular'].all():  ##different colors for circular crevasses and other crevasses
    #                 col = Col_Crevasses_Circ
    #             else:
    #                 continue ##We don't plot non circular crevasses anymre. Before we had: col = Col_Crevasses_Other
    #             ax4.plot(Df_plot['X'].values/1000, Df_plot['Y'].values/1000, color=col, linestyle='-', linewidth=2)
    #         ###Plot cavity contour
    #         ax4.plot(cavity_contour[:, 0]/1000,cavity_contour[:, 1]/1000,color=Col_Cavity,linestyle='-',linewidth=2.5)
    #         ###Make a horizontal color bar common to all subplots
    # cbar4bis = fig4bis.colorbar(CS4,
    #                     ax=axes4bis.ravel().tolist(),  # all subplots
    #                     orientation='horizontal',
    #                     ticks=levs_ticks4,
    #                     fraction=0.04,
    #                     pad=0.08)
    # cbar4bis.set_label(r'$\sigma_\mathrm{I}$ [kPa]', fontsize=20)
    # cbar4bis.ax.tick_params(labelsize=18)  # Adjust tick label size


    ###### START PLOT FIGURE 5 ######
    fig5.suptitle(r'Linear Visco-elastic $\nu = 0.45$', fontsize=23, weight='bold')
    print('starting Figure 5')
    for i, law in enumerate(['Subhour', 'Subday']):
        for j, tsp in enumerate(['a', 'b', 'c', 'd']):
            ##get proper ax
            ax5 = axes5[i, j]
            ###Get data for each subplots one by one
            df_cav5 = pd.DataFrame()
            df_nocav5 = pd.DataFrame()
            if i == 0:
                if j == 0:  ##Linear elastic
                    ax5.set_ylabel(r'Y [km]', fontsize=22)
                    df_nocav5 = Data_Simu_El[Data_Simu_El['Case'] == 'NoCavity']
                    df_cav5 = Data_Simu_El[Data_Simu_El['Case'] == 'EmptyCavity']
                    Title = r'Elastic $\nu = 0.3$'
                elif j == 1:  ##Linear visco elastic after 1s tsp =1s
                    df_nocav5 = Data_Simu_LVELCOMP_1s[(Data_Simu_LVELCOMP_1s['Timestep'] == 1) & (Data_Simu_LVELCOMP_1s['Case'] == 'NoCavity')]
                    df_cav5 = Data_Simu_LVELCOMP_1s[(Data_Simu_LVELCOMP_1s['Timestep'] == 1) & (Data_Simu_LVELCOMP_1s['Case'] == 'EmptyCavity')]
                    Title = r'@1s tsp = 1s'
                elif j == 2:  ##Linear visco elastic after 1mn tsp =1s
                    df_nocav5 = Data_Simu_LVELCOMP_1s[(Data_Simu_LVELCOMP_1s['Timestep'] == 60) & (Data_Simu_LVELCOMP_1s['Case'] == 'NoCavity')]
                    df_cav5 = Data_Simu_LVELCOMP_1s[(Data_Simu_LVELCOMP_1s['Timestep'] == 60) & (Data_Simu_LVELCOMP_1s['Case'] == 'EmptyCavity')]
                    Title = r'@1mn tsp = 1s'
                elif j == 3:  ##Linear visco elastic after 1h tsp=1mn
                    df_nocav5 = Data_Simu_LVELCOMP_5mn[(Data_Simu_LVELCOMP_5mn['Timestep'] == 12) & (Data_Simu_LVELCOMP_5mn['Case'] == 'NoCavity')]
                    df_cav5 = Data_Simu_LVELCOMP_5mn[(Data_Simu_LVELCOMP_5mn['Timestep'] == 12) & (Data_Simu_LVELCOMP_5mn['Case'] == 'EmptyCavity')]
                    Title = r'@1h tsp = 5mn'
            elif i == 1:
                ax5.set_xlabel(r'X [km]', fontsize=22)
                if j == 0:  ##Linear visco elastic after 6h tsp =5mn
                    ax5.set_ylabel(r'Y [km]', fontsize=22)
                    df_nocav5 = Data_Simu_LVELCOMP_5mn[(Data_Simu_LVELCOMP_5mn['Timestep'] == 6*60/5) & (Data_Simu_LVELCOMP_5mn['Case'] == 'NoCavity')]
                    df_cav5 = Data_Simu_LVELCOMP_5mn[(Data_Simu_LVELCOMP_5mn['Timestep'] == 6*60/5) & (Data_Simu_LVELCOMP_5mn['Case'] == 'EmptyCavity')]
                    Title = r'@6h tsp = 5mn'
                elif j == 1:  ##Linear visco elastic tsp = 12h tsp =5mn
                    df_nocav5 = Data_Simu_LVELCOMP_5mn[(Data_Simu_LVELCOMP_5mn['Timestep'] == 12*60/5) & (Data_Simu_LVELCOMP_5mn['Case'] == 'NoCavity')]
                    df_cav5 = Data_Simu_LVELCOMP_5mn[(Data_Simu_LVELCOMP_5mn['Timestep'] == 12*60/5) & (Data_Simu_LVELCOMP_5mn['Case'] == 'EmptyCavity')]
                    Title = r'@12h tsp = 5mn'
                elif j == 2:  ##Linear visco elastic tsp = 24h tsp =5mn
                    df_nocav5 = Data_Simu_LVELCOMP_5mn[(Data_Simu_LVELCOMP_5mn['Timestep'] == 24*60/5) & (Data_Simu_LVELCOMP_5mn['Case'] == 'NoCavity')]
                    df_cav5 = Data_Simu_LVELCOMP_5mn[(Data_Simu_LVELCOMP_5mn['Timestep'] == 24*60/5) & (Data_Simu_LVELCOMP_5mn['Case'] == 'EmptyCavity')]
                    Title = r'@24h tsp = 5mn'
                elif j == 3:  ##pure linear viscous incomp
                    df_nocav5 = Data_Simu_VL[Data_Simu_VL['Case'] == 'NoCavity']
                    df_cav5 = Data_Simu_VL[Data_Simu_VL['Case'] == 'EmptyCavity']
                    Title = r'pure viscous incomp.'
            ###Fill up the plots
            ax5.set_title(Title, fontsize=18, weight='bold')
            ax5.tick_params(labelsize=18)  # fontsize of the tick labels
            ax5.grid(True)
            ax5.grid(alpha=0.5)
            ###Force x and y axis to same scale
            ax5.set_aspect('equal', adjustable='box')
            ###SigmaI is max principal stress and at the surface always one of the principal stress is zero
            ### Therefore SigmaI is either positive or zero. Negative values of SigmaI are numerical artefacts
            df_nocav5.loc[df_nocav5["SigmaI"] < 0, "SigmaI"] = 0
            df_cav5.loc[df_cav5["SigmaI"] < 0, "SigmaI"] = 0
            ###Poject SigmaI on interpolation grid
            SigmaI_NoCav = Interpolate_field(df_nocav5, 'SigmaI', X, Y)
            SigmaI_Cav = Interpolate_field(df_cav5, 'SigmaI', X, Y)
            ###shading
            clevs5 = clevtight  ## cbar for shading
            # colorbar
            levs_ticks5 = levs_ticktight
            ###Fills up the map with colors for SigmaEq
            CS5 = ax5.contourf(X / 1000, Y / 1000, (SigmaI_Cav - SigmaI_NoCav) * 1000, clevs5, cmap=cmap,
                               extend='both')
            ax5.plot(xc / 1000, yc / 1000, color='k', linewidth=2)

            ###Below we remove colors that are outside of glacier contour
            clippath = mpltPath(np.c_[xc / 1000, yc / 1000])
            patch = PathPatch(clippath, facecolor='none')
            # ax = plt.gca()
            ax5.add_patch(patch)
            for c in [CS5]:
                c.set_clip_path(patch)
            ##plot crevasses as continuous line
            for crev_num in np.arange(Df_Crevasses['Crevasse Number'].min(),
                                      Df_Crevasses['Crevasse Number'].max() + 1):
                ###get proper points
                Df_plot = Df_Crevasses[Df_Crevasses['Crevasse Number'] == crev_num]
                if Df_plot['IsCircular'].all():  ##different colors for circular crevasses and other crevasses
                    col = Col_Crevasses_Circ
                else:
                    continue  ##We don't plot non circular crevasses anymre. Before we had: col = Col_Crevasses_Other
                ax5.plot(Df_plot['X'].values / 1000, Df_plot['Y'].values / 1000, color=col, linestyle='-',linewidth=2)
                ###Plot cavity contour
                ax5.plot(cavity_contour[:, 0] / 1000, cavity_contour[:, 1] / 1000, color=Col_Cavity, linestyle='-',linewidth=2.5)
            ###Make a horizontal color bar common to all subplots
    cbar5 = fig5.colorbar(CS5,
                          ax=axes5.ravel().tolist(),  # all subplots
                          orientation='horizontal',
                          ticks=levs_ticks5,
                          fraction=0.04,
                          pad=0.08)
    cbar5.set_label(r'$\Delta \sigma_\mathrm{I}$ [kPa]', fontsize=20)
    cbar5.ax.tick_params(labelsize=18)  # Adjust tick label size

    ###### START PLOT FIGURE 5 BIS######
    fig5bis.suptitle(r'Non Linear Visco-elastic $\nu = 0.45$', fontsize=23, weight='bold')
    print('starting Figure 5 bis')
    for i, law in enumerate(['Subhour', 'Subday']):
        for j, tsp in enumerate(['a', 'b', 'c', 'd']):
            ##get proper ax
            ax5 = axes5bis[i, j]
            ###Get data for each subplots one by one
            df_cav5 = pd.DataFrame()
            df_nocav5 = pd.DataFrame()
            if i == 0:
                if j == 0:  ##Linear elastic
                    ax5.set_ylabel(r'Y [km]', fontsize=22)
                    df_nocav5 = Data_Simu_El[Data_Simu_El['Case'] == 'NoCavity']
                    df_cav5 = Data_Simu_El[Data_Simu_El['Case'] == 'EmptyCavity']
                    Title = r'Elastic $\nu = 0.3$'
                elif j == 1:  ##Linear visco elastic after 1s tsp =1s
                    df_nocav5 = Data_Simu_NLVELCOMP_1s[(Data_Simu_NLVELCOMP_1s['Timestep'] == 1) & (Data_Simu_NLVELCOMP_1s['Case'] == 'NoCavity')]
                    df_cav5 = Data_Simu_NLVELCOMP_1s[(Data_Simu_NLVELCOMP_1s['Timestep'] == 1) & (Data_Simu_NLVELCOMP_1s['Case'] == 'EmptyCavity')]
                    Title = r'@1s tsp = 1s'
                elif j == 2:  ##Linear visco elastic after 1mn tsp =1s
                    df_nocav5 = Data_Simu_NLVELCOMP_1s[(Data_Simu_NLVELCOMP_1s['Timestep'] == 60) & (Data_Simu_NLVELCOMP_1s['Case'] == 'NoCavity')]
                    df_cav5 = Data_Simu_NLVELCOMP_1s[(Data_Simu_NLVELCOMP_1s['Timestep'] == 60) & (Data_Simu_NLVELCOMP_1s['Case'] == 'EmptyCavity')]
                    Title = r'@1mn tsp = 1s'
                elif j == 3:  ##Linear visco elastic after 1h tsp=1mn
                    df_nocav5 = Data_Simu_NLVELCOMP_5mn[(Data_Simu_NLVELCOMP_5mn['Timestep'] == 12) & (Data_Simu_NLVELCOMP_5mn['Case'] == 'NoCavity')]
                    df_cav5 = Data_Simu_NLVELCOMP_5mn[(Data_Simu_NLVELCOMP_5mn['Timestep'] == 12) & (Data_Simu_NLVELCOMP_5mn['Case'] == 'EmptyCavity')]
                    Title = r'@1h tsp = 5mn'
            elif i == 1:
                ax5.set_xlabel(r'X [km]', fontsize=22)
                if j == 0:  ##Linear visco elastic tsp = 6h
                    ax5.set_ylabel(r'Y [km]', fontsize=22)
                    df_nocav5 = Data_Simu_NLVELCOMP_5mn[(Data_Simu_NLVELCOMP_5mn['Timestep'] == 6 * 60 / 5) & (Data_Simu_NLVELCOMP_5mn['Case'] == 'NoCavity')]
                    df_cav5 = Data_Simu_NLVELCOMP_5mn[(Data_Simu_NLVELCOMP_5mn['Timestep'] == 6 * 60 / 5) & (Data_Simu_NLVELCOMP_5mn['Case'] == 'EmptyCavity')]
                    Title = r'@6h tsp = 5mn'
                elif j == 1:  ##Linear visco elastic tsp = 12h tsp =5mn
                    df_nocav5 = Data_Simu_NLVELCOMP_5mn[(Data_Simu_NLVELCOMP_5mn['Timestep'] == 12 * 60 / 5) & (Data_Simu_NLVELCOMP_5mn['Case'] == 'NoCavity')]
                    df_cav5 = Data_Simu_NLVELCOMP_5mn[(Data_Simu_NLVELCOMP_5mn['Timestep'] == 12 * 60 / 5) & (Data_Simu_NLVELCOMP_5mn['Case'] == 'EmptyCavity')]
                    Title = r'@12h tsp = 5mn'
                elif j == 2:  ##Linear visco elastic tsp = 24h tsp =5mn
                    df_nocav5 = Data_Simu_NLVELCOMP_5mn[(Data_Simu_NLVELCOMP_5mn['Timestep'] == 24 * 60 / 5) & (Data_Simu_NLVELCOMP_5mn['Case'] == 'NoCavity')]
                    df_cav5 = Data_Simu_NLVELCOMP_5mn[(Data_Simu_NLVELCOMP_5mn['Timestep'] == 24 * 60 / 5) & (Data_Simu_NLVELCOMP_5mn['Case'] == 'EmptyCavity')]
                    Title = r'@24h tsp = 5mn'
                elif j == 3:  ##pure Glen
                    df_nocav5 = Data_Simu_VNL[Data_Simu_VNL['Case'] == 'NoCavity']
                    df_cav5 = Data_Simu_VNL[Data_Simu_VNL['Case'] == 'EmptyCavity']
                    Title = 'Glen-Nye incomp.'
            ###Fill up the plots
            ax5.set_title(Title, fontsize=18, weight='bold')
            ax5.tick_params(labelsize=18)  # fontsize of the tick labels
            ax5.grid(True)
            ax5.grid(alpha=0.5)
            ###Force x and y axis to same scale
            ax5.set_aspect('equal', adjustable='box')
            ###SigmaI is max principal stress and at the surface always one of the principal stress is zero
            ### Therefore SigmaI is either positive or zero. Negative values of SigmaI are numerical artefacts
            df_nocav5.loc[df_nocav5["SigmaI"] < 0, "SigmaI"] = 0
            df_cav5.loc[df_cav5["SigmaI"] < 0, "SigmaI"] = 0
            ###Poject SigmaI on interpolation grid
            SigmaI_NoCav = Interpolate_field(df_nocav5, 'SigmaI', X, Y)
            SigmaI_Cav = Interpolate_field(df_cav5, 'SigmaI', X, Y)
            ###shading
            clevs5 = clevwide ## cbar for shading
            # colorbar
            levs_ticks5 = levs_tickwide
            ###Fills up the map with colors for SigmaEq
            CS5 = ax5.contourf(X / 1000, Y / 1000, (SigmaI_Cav - SigmaI_NoCav) * 1000, clevs5, cmap=cmap,
                               extend='both')
            ax5.plot(xc / 1000, yc / 1000, color='k', linewidth=2)

            ###Below we remove colors that are outside of glacier contour
            clippath = mpltPath(np.c_[xc / 1000, yc / 1000])
            patch = PathPatch(clippath, facecolor='none')
            # ax = plt.gca()
            ax5.add_patch(patch)
            for c in [CS5]:
                c.set_clip_path(patch)
            ##plot crevasses as continuous line
            for crev_num in np.arange(Df_Crevasses['Crevasse Number'].min(),
                                      Df_Crevasses['Crevasse Number'].max() + 1):
                ###get proper points
                Df_plot = Df_Crevasses[Df_Crevasses['Crevasse Number'] == crev_num]
                if Df_plot['IsCircular'].all():  ##different colors for circular crevasses and other crevasses
                    col = Col_Crevasses_Circ
                else:
                    continue  ##We don't plot non circular crevasses anymre. Before we had: col = Col_Crevasses_Other
                ax5.plot(Df_plot['X'].values / 1000, Df_plot['Y'].values / 1000, color=col, linestyle='-',
                         linewidth=2)
            ###Plot cavity contour
            ax5.plot(cavity_contour[:, 0] / 1000, cavity_contour[:, 1] / 1000, color=Col_Cavity, linestyle='-',linewidth=2.5)
            ###Make a horizontal color bar common to all subplots
    cbar5bis = fig5bis.colorbar(CS5,
                          ax=axes5bis.ravel().tolist(),  # all subplots
                          orientation='horizontal',
                          ticks=levs_ticks5,
                          fraction=0.04,
                          pad=0.08)
    cbar5bis.set_label(r'$\Delta \sigma_\mathrm{I}$ [kPa]', fontsize=20)
    cbar5bis.ax.tick_params(labelsize=18)  # Adjust tick label size

    # ###### START PLOT FIGURE 6 ######
    # fig6.suptitle('Sensitivity to timestep (Empty Cavity @5mn)', fontsize=23, weight='bold')
    # print('starting Figure 6')
    # for i, law in enumerate(['LVELCOMP', 'LVELINCOMP', 'NLVELCOMP', 'NLVELINCOMP']):
    #     df_1s = pd.DataFrame()
    #     df_1mn = pd.DataFrame()
    #     if law == 'LVELCOMP':
    #         df_1s = Data_Simu_LVELCOMP_1s[(Data_Simu_LVELCOMP_1s['Timestep'] == 5*60) & (Data_Simu_LVELCOMP_1s['Case'] == 'EmptyCavity')]
    #         df_1mn = Data_Simu_LVELCOMP_1mn[(Data_Simu_LVELCOMP_1mn['Timestep'] == 5) & (Data_Simu_LVELCOMP_1mn['Case'] == 'EmptyCavity')]
    #         linetitle = r'Li. $\nu =0.3$'
    #     elif law == 'LVELINCOMP':
    #         df_1s = Data_Simu_LVELINCOMP_1s[(Data_Simu_LVELINCOMP_1s['Timestep'] == 5*60) & (Data_Simu_LVELINCOMP_1s['Case'] == 'EmptyCavity')]
    #         df_1mn = Data_Simu_LVELINCOMP_1mn[(Data_Simu_LVELINCOMP_1mn['Timestep'] == 5) & (Data_Simu_LVELINCOMP_1mn['Case'] == 'EmptyCavity')]
    #         linetitle = r'Li. incomp.'
    #     elif law == 'NLVELCOMP':
    #         df_1s = Data_Simu_NLVELCOMP_1s[(Data_Simu_NLVELCOMP_1s['Timestep'] == 5*60) & (Data_Simu_NLVELCOMP_1s['Case'] == 'EmptyCavity')]
    #         df_1mn = Data_Simu_NLVELCOMP_1mn[(Data_Simu_NLVELCOMP_1mn['Timestep'] == 5) & (Data_Simu_NLVELCOMP_1mn['Case'] == 'EmptyCavity')]
    #         linetitle = r'NL $\nu =0.3$'
    #     elif law == 'NLVELINCOMP':
    #         df_1s = Data_Simu_NLVELINCOMP_1s[(Data_Simu_NLVELINCOMP_1s['Timestep'] == 5*60) & (Data_Simu_NLVELINCOMP_1s['Case'] == 'EmptyCavity')]
    #         df_1mn= Data_Simu_NLVELINCOMP_1mn[(Data_Simu_NLVELINCOMP_1mn['Timestep'] == 5) & (Data_Simu_NLVELINCOMP_1mn['Case'] == 'EmptyCavity')]
    #         linetitle = r'NL incomp.'
    #     ###Poject SigmaI on interpolation grid
    #     SigmaI_tsp1s = Interpolate_field(df_1s, 'SigmaI', X, Y)
    #     SigmaI_tsp1mn = Interpolate_field(df_1mn, 'SigmaI', X, Y)
    #     SigmaI_diff = SigmaI_tsp1mn - SigmaI_tsp1s
    #     for j,tsp in enumerate(['1s','1mn','diff']):
    #         ax6 = axes6[i, j]
    #         if i == 3:
    #             ax6.set_xlabel(r'X [km]', fontsize=22)
    #         if j == 0:
    #             ax6.set_ylabel(r'Y [km]', fontsize=22)
    #         ax6.tick_params(labelsize=18)  # fontsize of the tick labels
    #         ax6.grid(alpha=0.5)
    #         ###Force x and y axis to same scale
    #         ax6.set_aspect('equal', adjustable='box')
    #         ###shading
    #         if j<2:
    #             clevs6 = clevwide  ## cbar for shading
    #             # colorbar
    #             levs_ticks6 = levs_tickwide
    #         else:
    #             if AbsoluteDifference:
    #                 clevs6 = np.arange(-0.05001, 0.05001, 0.01)  ## cbar for shading
    #                 # colorbar
    #                 levs_ticks6diff = np.arange(-0.050001, 0.050001, 0.05)
    #             else:
    #                 clevs6 = np.arange(-2.0001, 2.0001, 0.02)  ## cbar for shading
    #                 # colorbar
    #                 levs_ticks6diff = np.arange(-2.0001, 2.0001, 1.0)
    #         ###Fills up the map with colors for SigmaEq
    #         if tsp == '1s':
    #             CS6 = ax6.contourf(X / 1000, Y / 1000, (SigmaI_tsp1s) * 1000, clevs6, cmap=cmap, extend='both')
    #             Title = r'tsp = 1s'
    #         elif tsp == '1mn':
    #             CS6 = ax6.contourf(X / 1000, Y / 1000, (SigmaI_tsp1mn) * 1000, clevs6, cmap=cmap, extend='both')
    #             Title = r'tsp = 1mn'
    #         elif tsp == 'diff':
    #             if AbsoluteDifference:
    #                 CS6_diff = ax6.contourf(X / 1000, Y / 1000, SigmaI_diff, clevs6, cmap=cmapdiff, extend='both')
    #             else:
    #                 CS6_diff = ax6.contourf(X / 1000, Y / 1000, (SigmaI_diff) / SigmaI_tsp1s, clevs6, cmap=cmapdiff, extend='both')
    #             Title = r'Diff: 1mn - 1s'
    #         if i==0:
    #             ax6.set_title(Title, fontsize=16)
    #         if j ==0:
    #             ax6.text(-0.6, 0.5, linetitle, transform=ax6.transAxes, rotation=90, ha='center', va='center', fontsize=14, weight='bold')
    #         ax6.plot(xc / 1000, yc / 1000, color='k', linewidth=2)
    #         ###Below we remove colors that are outside of glacier contour
    #         clippath = mpltPath(np.c_[xc / 1000, yc / 1000])
    #         patch = PathPatch(clippath, facecolor='none')
    #         # ax = plt.gca()
    #         ax6.add_patch(patch)
    #         if j<2:
    #             for c in [CS6]:
    #                 c.set_clip_path(patch)
    #         else:
    #             for c in [CS6_diff]:
    #                 c.set_clip_path(patch)
    #         ##plot crevasses as continuous line
    #         for crev_num in np.arange(Df_Crevasses['Crevasse Number'].min(), Df_Crevasses['Crevasse Number'].max() + 1):
    #             ###get proper points
    #             Df_plot = Df_Crevasses[Df_Crevasses['Crevasse Number'] == crev_num]
    #             if Df_plot['IsCircular'].all():  ##different colors for circular crevasses and other crevasses
    #                 col = Col_Crevasses_Circ
    #             else:
    #                 continue  ##We don't plot non circular crevasses anymre. Before we had: col = Col_Crevasses_Other
    #             ax6.plot(Df_plot['X'].values / 1000, Df_plot['Y'].values / 1000, color=col, linestyle='-', linewidth=2)
    #         ###Plot cavity contour
    #         ax6.plot(cavity_contour[:, 0] / 1000, cavity_contour[:, 1] / 1000, color=Col_Cavity, linestyle='-',linewidth=2.5)
    #
    # # After the loop (once everything is plotted):
    # cbar1 = fig6.colorbar(CS6, ax=axes6[:, :2].ravel().tolist(), orientation='horizontal', ticks=levs_ticks6, fraction=0.04, pad=0.08 )
    # cbar1.set_label(r'$\sigma_\mathrm{I}$ [kPa]', fontsize=20)
    # cbar1.ax.tick_params(labelsize=18)
    # cbar2 = fig6.colorbar( CS6_diff, ax=axes6[:, 2].ravel().tolist(), orientation='horizontal', ticks=levs_ticks6diff, fraction=0.08, pad=0.08 )
    # if AbsoluteDifference:
    #     cbar2.set_label(r'$\Delta\sigma_\mathrm{I}$ [kPa] ', fontsize=20)
    # else:
    #     cbar2.set_label(r'$\Delta\sigma_\mathrm{I}/ \sigma_\mathrm{I}$ ', fontsize=20)
    # cbar2.ax.tick_params(labelsize=18)


    # ###### START PLOT FIGURE 6 BIS ######
    # fig6bis.suptitle('Sensitivity to timestep (No Cavity @10mn)', fontsize=23, weight='bold')
    # print('starting Figure 6 bis')
    # for i, law in enumerate(['LVELCOMP', 'LVELINCOMP', 'NLVELCOMP', 'NLVELINCOMP']):
    #     df_1mn = pd.DataFrame()
    #     df_5mn = pd.DataFrame()
    #     if law == 'LVELCOMP':
    #         df_1mn = Data_Simu_LVELCOMP_1mn[(Data_Simu_LVELCOMP_1mn['Timestep'] == 10) & (Data_Simu_LVELCOMP_1mn['Case'] == 'EmptyCavity')]
    #         df_5mn = Data_Simu_LVELCOMP_5mn[(Data_Simu_LVELCOMP_5mn['Timestep'] == 2) & (Data_Simu_LVELCOMP_5mn['Case'] == 'EmptyCavity')]
    #         linetitle = r'Li. $\nu =0.3$'
    #     elif law == 'LVELINCOMP':
    #         df_1mn = Data_Simu_LVELINCOMP_1mn[(Data_Simu_LVELINCOMP_1mn['Timestep'] == 10) & (Data_Simu_LVELINCOMP_1mn['Case'] == 'EmptyCavity')]
    #         df_5mn = Data_Simu_LVELINCOMP_5mn[(Data_Simu_LVELINCOMP_5mn['Timestep'] == 2) & (Data_Simu_LVELINCOMP_5mn['Case'] == 'EmptyCavity')]
    #         linetitle = r'Li. incomp.'
    #     elif law == 'NLVELCOMP':
    #         df_1mn = Data_Simu_NLVELCOMP_1mn[(Data_Simu_NLVELCOMP_1mn['Timestep'] == 10) & (Data_Simu_NLVELCOMP_1mn['Case'] == 'EmptyCavity')]
    #         df_5mn = Data_Simu_NLVELCOMP_5mn[(Data_Simu_NLVELCOMP_5mn['Timestep'] == 2) & (Data_Simu_NLVELCOMP_5mn['Case'] == 'EmptyCavity')]
    #         linetitle = r'NL $\nu =0.3$'
    #     elif law == 'NLVELINCOMP':
    #         df_1mn = Data_Simu_NLVELINCOMP_1mn[(Data_Simu_NLVELINCOMP_1mn['Timestep'] == 10) & (Data_Simu_NLVELINCOMP_1mn['Case'] == 'EmptyCavity')]
    #         df_5mn = Data_Simu_NLVELINCOMP_5mn[(Data_Simu_NLVELINCOMP_5mn['Timestep'] == 2) & (Data_Simu_NLVELINCOMP_5mn['Case'] == 'EmptyCavity')]
    #         linetitle = r'NL incomp.'
    #     ###Poject SigmaI on interpolation grid
    #     SigmaI_tsp1mn = Interpolate_field(df_1mn, 'SigmaI', X, Y)
    #     SigmaI_tsp5mn = Interpolate_field(df_5mn, 'SigmaI', X, Y)
    #     SigmaI_diff = SigmaI_tsp5mn - SigmaI_tsp1mn
    #     for j,tsp in enumerate(['1mn','5mn','diff']):
    #         ax6 = axes6bis[i, j]
    #         if i == 3:
    #             ax6.set_xlabel(r'X [km]', fontsize=22)
    #         if j == 0:
    #             ax6.set_ylabel(r'Y [km]', fontsize=22)
    #         ax6.tick_params(labelsize=18)  # fontsize of the tick labels
    #         ax6.grid(alpha=0.5)
    #         ###Force x and y axis to same scale
    #         ax6.set_aspect('equal', adjustable='box')
    #         ###shading
    #         if j<2:
    #             clevs6 = clevwide  ## cbar for shading
    #             # colorbar
    #             levs_ticks6 = levs_tickwide
    #         else:
    #             if AbsoluteDifference:
    #                 clevs6 = np.arange(-0.0501, 0.0501, 0.01)  ## cbar for shading
    #                 # colorbar
    #                 levs_ticks6diff = np.arange(-0.050001, 0.050001, 0.05)
    #             else:
    #                 clevs6 = np.arange(-2.0001, 2.0001, 0.02)  ## cbar for shading
    #                 # colorbar
    #                 levs_ticks6diff = np.arange(-2.0001, 2.0001, 1.0)
    #         ###Fills up the map with colors for SigmaEq
    #         if tsp == '1mn':
    #             CS6 = ax6.contourf(X / 1000, Y / 1000, (SigmaI_tsp1mn) * 1000, clevs6, cmap=cmap, extend='both')
    #             Title = r'tsp = 1mn'
    #         elif tsp == '5mn':
    #             CS6 = ax6.contourf(X / 1000, Y / 1000, (SigmaI_tsp5mn) * 1000, clevs6, cmap=cmap, extend='both')
    #             Title = r'tsp = 5mn'
    #         elif tsp == 'diff':
    #             if AbsoluteDifference:
    #                 CS6_diff = ax6.contourf(X / 1000, Y / 1000, SigmaI_diff, clevs6, cmap=cmapdiff, extend='both')
    #                 Title = r'Abs. Diff: 5mn - 1mn'
    #             else:
    #                 CS6_diff = ax6.contourf(X / 1000, Y / 1000, (SigmaI_diff) / SigmaI_tsp1mn, clevs6, cmap=cmapdiff, extend='both')
    #                 Title = r'Rel. Diff: 5mn - 1mn'
    #         if i==0:
    #             ax6.set_title(Title, fontsize=16)
    #         if j ==0:
    #             ax6.text(-0.6, 0.5, linetitle, transform=ax6.transAxes, rotation=90, ha='center', va='center', fontsize=14, weight='bold')
    #         ax6.plot(xc / 1000, yc / 1000, color='k', linewidth=2)
    #         ###Below we remove colors that are outside of glacier contour
    #         clippath = mpltPath(np.c_[xc / 1000, yc / 1000])
    #         patch = PathPatch(clippath, facecolor='none')
    #         # ax = plt.gca()
    #         ax6.add_patch(patch)
    #         if j<2:
    #             for c in [CS6]:
    #                 c.set_clip_path(patch)
    #         else:
    #             for c in [CS6_diff]:
    #                 c.set_clip_path(patch)
    #         ##plot crevasses as continuous line
    #         for crev_num in np.arange(Df_Crevasses['Crevasse Number'].min(), Df_Crevasses['Crevasse Number'].max() + 1):
    #             ###get proper points
    #             Df_plot = Df_Crevasses[Df_Crevasses['Crevasse Number'] == crev_num]
    #             if Df_plot['IsCircular'].all():  ##different colors for circular crevasses and other crevasses
    #                 col = Col_Crevasses_Circ
    #             else:
    #                 continue  ##We don't plot non circular crevasses anymre. Before we had: col = Col_Crevasses_Other
    #             ax6.plot(Df_plot['X'].values / 1000, Df_plot['Y'].values / 1000, color=col, linestyle='-', linewidth=2)
    #         ###Plot cavity contour
    #         ax6.plot(cavity_contour[:, 0] / 1000, cavity_contour[:, 1] / 1000, color=Col_Cavity, linestyle='-',linewidth=2.5)
    #
    # # After the loop (once everything is plotted):
    # cbar1bis = fig6bis.colorbar(CS6, ax=axes6bis[:, :2].ravel().tolist(), orientation='horizontal', ticks=levs_ticks6, fraction=0.04, pad=0.08 )
    # cbar1bis.set_label(r'$\sigma_\mathrm{I}$ [kPa]', fontsize=20)
    # cbar1bis.ax.tick_params(labelsize=18)
    # cbar2bis = fig6bis.colorbar( CS6_diff, ax=axes6bis[:, 2].ravel().tolist(), orientation='horizontal', ticks=levs_ticks6diff, fraction=0.08, pad=0.08 )
    # if AbsoluteDifference:
    #     cbar2bis.set_label(r'$\Delta\sigma_\mathrm{I}$ [kPa] ', fontsize=20)
    # else:
    #     cbar2bis.set_label(r'$\Delta\sigma_\mathrm{I}/ \sigma_\mathrm{I}$ ', fontsize=20)
    # cbar2bis.ax.tick_params(labelsize=18)
    #
    # ###### START PLOT FIGURE 7 ######
    # fig7.suptitle('Sensitivity to timestep (No Cavity @5h)', fontsize=23, weight='bold')
    # print('starting Figure 7')
    # for i, law in enumerate(['LVELCOMP', 'LVELINCOMP', 'NLVELCOMP', 'NLVELINCOMP']):
    #     df_1mn = pd.DataFrame()
    #     df_5mn = pd.DataFrame()
    #     if law == 'LVELCOMP':
    #         df_1mn = Data_Simu_LVELCOMP_1mn[(Data_Simu_LVELCOMP_1mn['Timestep'] == 5*60) & (Data_Simu_LVELCOMP_1mn['Case'] == 'EmptyCavity')]
    #         df_5mn = Data_Simu_LVELCOMP_5mn[(Data_Simu_LVELCOMP_5mn['Timestep'] == 5*12) & (Data_Simu_LVELCOMP_5mn['Case'] == 'EmptyCavity')]
    #         linetitle = r'Li. $\nu =0.3$'
    #     elif law == 'LVELINCOMP':
    #         df_1mn = Data_Simu_LVELINCOMP_1mn[(Data_Simu_LVELINCOMP_1mn['Timestep'] == 5*60) & (Data_Simu_LVELINCOMP_1mn['Case'] == 'EmptyCavity')]
    #         df_5mn = Data_Simu_LVELINCOMP_5mn[(Data_Simu_LVELINCOMP_5mn['Timestep'] == 5*12) & (Data_Simu_LVELINCOMP_5mn['Case'] == 'EmptyCavity')]
    #         linetitle = r'Li. incomp.'
    #     elif law == 'NLVELCOMP':
    #         df_1mn = Data_Simu_NLVELCOMP_1mn[(Data_Simu_NLVELCOMP_1mn['Timestep'] == 5*60) & (Data_Simu_NLVELCOMP_1mn['Case'] == 'EmptyCavity')]
    #         df_5mn = Data_Simu_NLVELCOMP_5mn[(Data_Simu_NLVELCOMP_5mn['Timestep'] == 5*12) & (Data_Simu_NLVELCOMP_5mn['Case'] == 'EmptyCavity')]
    #         linetitle = r'NL $\nu =0.3$'
    #     elif law == 'NLVELINCOMP':
    #         df_1mn = Data_Simu_NLVELINCOMP_1mn[(Data_Simu_NLVELINCOMP_1mn['Timestep'] == 5*60) & (Data_Simu_NLVELINCOMP_1mn['Case'] == 'EmptyCavity')]
    #         df_5mn = Data_Simu_NLVELINCOMP_5mn[(Data_Simu_NLVELINCOMP_5mn['Timestep'] == 5*12) & (Data_Simu_NLVELINCOMP_5mn['Case'] == 'EmptyCavity')]
    #         linetitle = r'NL incomp.'
    #     ###Poject SigmaI on interpolation grid
    #     SigmaI_tsp5mn = Interpolate_field(df_5mn, 'SigmaI', X, Y)
    #     SigmaI_tsp1mn = Interpolate_field(df_1mn, 'SigmaI', X, Y)
    #     SigmaI_diff = SigmaI_tsp5mn - SigmaI_tsp1mn
    #     for j,tsp in enumerate(['1mn','5mn','diff']):
    #         ax7 = axes7[i, j]
    #         if i == 3:
    #             ax7.set_xlabel(r'X [km]', fontsize=22)
    #         if j == 0:
    #             ax7.set_ylabel(r'Y [km]', fontsize=22)
    #         ax7.tick_params(labelsize=18)  # fontsize of the tick labels
    #         ax7.grid(alpha=0.5)
    #         ###Force x and y axis to same scale
    #         ax7.set_aspect('equal', adjustable='box')
    #         ###shading
    #         if j<2:
    #             clevs7 = clevwide  ## cbar for shading
    #             # colorbar
    #             levs_ticks7 = levs_tickwide
    #         else:
    #             if AbsoluteDifference:
    #                 clevs7 = np.arange(-0.0501, 0.0501, 0.01)  ## cbar for shading
    #                 # colorbar
    #                 levs_ticks7diff = np.arange(-0.050001, 0.050001, 0.05)
    #             else:
    #                 clevs7 = np.arange(-2.0, 2.0, 0.02)  ## cbar for shading
    #                 # colorbar
    #                 levs_ticks7diff = np.arange(-2.0001, 2.0001, 1.0)
    #         ###Fills up the map with colors for SigmaEq
    #         if tsp == '1mn':
    #             CS7 = ax7.contourf(X / 1000, Y / 1000, (SigmaI_tsp1mn) * 1000, clevs7, cmap=cmap, extend='both')
    #             Title = r'tsp = 1mn'
    #         elif tsp == '5mn':
    #             CS7 = ax7.contourf(X / 1000, Y / 1000, (SigmaI_tsp5mn) * 1000, clevs7, cmap=cmap, extend='both')
    #             Title = r'tsp = 5mn'
    #         elif tsp == 'diff':
    #             if AbsoluteDifference:
    #                 CS7_diff = ax7.contourf(X / 1000, Y / 1000, (SigmaI_diff), clevs7, cmap=cmapdiff, extend='both')
    #                 Title = r'Abs. Diff: 5mn - 1mn'
    #             else:
    #                 CS7_diff = ax7.contourf(X / 1000, Y / 1000, (SigmaI_diff) /(SigmaI_tsp1mn), clevs7, cmap=cmapdiff, extend='both')
    #                 Title = r'Rel. Diff: 5mn - 1mn'
    #         if i==0:
    #             ax7.set_title(Title, fontsize=16)
    #         if j ==0:
    #             ax7.text(-0.6, 0.5, linetitle, transform=ax7.transAxes, rotation=90, ha='center', va='center', fontsize=14, weight='bold')
    #         ax7.plot(xc / 1000, yc / 1000, color='k', linewidth=2)
    #         ###Below we remove colors that are outside of glacier contour
    #         clippath = mpltPath(np.c_[xc / 1000, yc / 1000])
    #         patch = PathPatch(clippath, facecolor='none')
    #         # ax = plt.gca()
    #         ax7.add_patch(patch)
    #         if j<2:
    #             for c in [CS7]:
    #                 c.set_clip_path(patch)
    #         else:
    #             for c in [CS7_diff]:
    #                 c.set_clip_path(patch)
    #         ##plot crevasses as continuous line
    #         for crev_num in np.arange(Df_Crevasses['Crevasse Number'].min(), Df_Crevasses['Crevasse Number'].max() + 1):
    #             ###get proper points
    #             Df_plot = Df_Crevasses[Df_Crevasses['Crevasse Number'] == crev_num]
    #             if Df_plot['IsCircular'].all():  ##different colors for circular crevasses and other crevasses
    #                 col = Col_Crevasses_Circ
    #             else:
    #                 continue  ##We don't plot non circular crevasses anymre. Before we had: col = Col_Crevasses_Other
    #             ax7.plot(Df_plot['X'].values / 1000, Df_plot['Y'].values / 1000, color=col, linestyle='-', linewidth=2)
    #         ###Plot cavity contour
    #         ax7.plot(cavity_contour[:, 0] / 1000, cavity_contour[:, 1] / 1000, color=Col_Cavity, linestyle='-',linewidth=2.5)
    #
    # # After the loop (once everything is plotted):
    # cbar1 = fig7.colorbar(CS7, ax=axes7[:, :2].ravel().tolist(), orientation='horizontal', ticks=levs_ticks7, fraction=0.04, pad=0.08 )
    # cbar1.set_label(r'$\sigma_\mathrm{I}$ [kPa]', fontsize=20)
    # cbar1.ax.tick_params(labelsize=18)
    # cbar2 = fig7.colorbar( CS7_diff, ax=axes7[:, 2].ravel().tolist(), orientation='horizontal', ticks=levs_ticks7diff, fraction=0.08, pad=0.08 )
    # if AbsoluteDifference:
    #     cbar2.set_label(r'$\Delta\sigma_\mathrm{I}$ [kPa]', fontsize=20)
    # else:
    #     cbar2.set_label(r'$\Delta\sigma_\mathrm{I}/ \sigma_\mathrm{I}$ ', fontsize=20)
    # cbar2.ax.tick_params(labelsize=18)
    #
    # #
    #

    #
    # ###### START PLOT FIGURE 8 ######
    # fig8.suptitle(r'pure Viscous Linear transient', fontsize=23, weight='bold')
    # print('starting Figure 8')
    # for i, law in enumerate(['Subhour', 'Subday']):
    #     for j, tsp in enumerate(['a', 'b', 'c', 'd']):
    #         ##get proper ax
    #         ax8 = axes8[i, j]
    #         ###Get data for each subplots one by one
    #         df_cav8 = pd.DataFrame()
    #         df_nocav8 = pd.DataFrame()
    #         if i == 0:
    #             if j == 0:  ##Linear elastic
    #                 ax8.set_ylabel(r'Y [km]', fontsize=22)
    #                 df_nocav8 = Data_Simu_El[Data_Simu_El['Case'] == 'NoCavity']
    #                 df_cav8 = Data_Simu_El[Data_Simu_El['Case'] == 'EmptyCavity']
    #                 Title = r'Elastic $\nu = 0.3$'
    #             elif j == 1:  ##Linear viscous after 1s tsp =1s
    #                 df_nocav8 = Data_Simu_VLTRANS_1s[(Data_Simu_VLTRANS_1s['Timestep'] == 1) & (Data_Simu_VLTRANS_1s['Case'] == 'NoCavity')]
    #                 df_cav8 = Data_Simu_VLTRANS_1s[(Data_Simu_VLTRANS_1s['Timestep'] == 1) & (Data_Simu_VLTRANS_1s['Case'] == 'EmptyCavity')]
    #                 Title = r'@1s tsp = 1s'
    #             elif j == 2:  ##Linear viscous after 1mn tsp =1s
    #                 df_nocav8 = Data_Simu_VLTRANS_1s[(Data_Simu_VLTRANS_1s['Timestep'] == 60) & (Data_Simu_VLTRANS_1s['Case'] == 'NoCavity')]
    #                 df_cav8 = Data_Simu_VLTRANS_1s[(Data_Simu_VLTRANS_1s['Timestep'] == 60) & (Data_Simu_VLTRANS_1s['Case'] == 'EmptyCavity')]
    #                 Title = r'@1mn tsp = 1s'
    #             elif j == 3:  ##Linear viscous after 1h tsp=1mn
    #                 df_nocav8 = Data_Simu_VLTRANS_5mn[(Data_Simu_VLTRANS_5mn['Timestep'] == 12) & (Data_Simu_VLTRANS_5mn['Case'] == 'NoCavity')]
    #                 df_cav8 = Data_Simu_VLTRANS_5mn[(Data_Simu_VLTRANS_5mn['Timestep'] == 12) & (Data_Simu_VLTRANS_5mn['Case'] == 'EmptyCavity')]
    #                 Title = r'@1h tsp = 5mn'
    #         elif i == 1:
    #             ax8.set_xlabel(r'X [km]', fontsize=22)
    #             if j == 0:  ##Linear viscous after 6h tsp =5mn
    #                 ax8.set_ylabel(r'Y [km]', fontsize=22)
    #                 df_nocav8 = Data_Simu_VLTRANS_5mn[(Data_Simu_VLTRANS_5mn['Timestep'] == 6*60/5) & (Data_Simu_VLTRANS_5mn['Case'] == 'NoCavity')]
    #                 df_cav8 = Data_Simu_VLTRANS_5mn[(Data_Simu_VLTRANS_5mn['Timestep'] == 6*60/5) & (Data_Simu_VLTRANS_5mn['Case'] == 'EmptyCavity')]
    #                 Title = r'@6h tsp = 5mn'
    #             elif j == 1:  ##Linear viscous tsp = 12h tsp =5mn
    #                 df_nocav8 = Data_Simu_VLTRANS_5mn[(Data_Simu_VLTRANS_5mn['Timestep'] == 12*60/5) & (Data_Simu_VLTRANS_5mn['Case'] == 'NoCavity')]
    #                 df_cav8 = Data_Simu_VLTRANS_5mn[(Data_Simu_VLTRANS_5mn['Timestep'] == 12*60/5) & (Data_Simu_VLTRANS_5mn['Case'] == 'EmptyCavity')]
    #                 Title = r'@12h tsp = 5mn'
    #             elif j == 2:  ##Linear viscous tsp = 24h tsp =5mn
    #                 df_nocav8 = Data_Simu_VLTRANS_5mn[(Data_Simu_VLTRANS_5mn['Timestep'] == 24*60/5) & (Data_Simu_VLTRANS_5mn['Case'] == 'NoCavity')]
    #                 df_cav8 = Data_Simu_VLTRANS_5mn[(Data_Simu_VLTRANS_5mn['Timestep'] == 24*60/5) & (Data_Simu_VLTRANS_5mn['Case'] == 'EmptyCavity')]
    #                 Title = r'@24h tsp = 5mn'
    #             elif j == 3:  ##pure linear viscous incomp
    #                 df_nocav8 = Data_Simu_VL[Data_Simu_VL['Case'] == 'NoCavity']
    #                 df_cav8 = Data_Simu_VL[Data_Simu_VL['Case'] == 'EmptyCavity']
    #                 Title = r'pure viscous steady'
    #         ###Fill up the plots
    #         ax8.set_title(Title, fontsize=18, weight='bold')
    #         ax8.tick_params(labelsize=18)  # fontsize of the tick labels
    #         ax8.grid(True)
    #         ax8.grid(alpha=0.5)
    #         ###Force x and y axis to same scale
    #         ax8.set_aspect('equal', adjustable='box')
    #         ###Poject SigmaI on interpolation grid
    #         SigmaI_NoCav = Interpolate_field(df_nocav8, 'SigmaI', X, Y)
    #         SigmaI_Cav = Interpolate_field(df_cav8, 'SigmaI', X, Y)
    #         ###shading
    #         clevs8 = clevtight  ## cbar for shading
    #         # colorbar
    #         levs_ticks8 = levs_ticktight
    #         ###Fills up the map with colors for SigmaEq
    #         CS8 = ax8.contourf(X / 1000, Y / 1000, (SigmaI_Cav - SigmaI_NoCav) * 1000, clevs8, cmap=cmap, extend='both')
    #         ax8.plot(xc / 1000, yc / 1000, color='k', linewidth=2)
    #         ###Below we remove colors that are outside of glacier contour
    #         clippath = mpltPath(np.c_[xc / 1000, yc / 1000])
    #         patch = PathPatch(clippath, facecolor='none')
    #         # ax = plt.gca()
    #         ax8.add_patch(patch)
    #         for c in [CS8]:
    #             c.set_clip_path(patch)
    #         ##plot crevasses as continuous line
    #         for crev_num in np.arange(Df_Crevasses['Crevasse Number'].min(),
    #                                   Df_Crevasses['Crevasse Number'].max() + 1):
    #             ###get proper points
    #             Df_plot = Df_Crevasses[Df_Crevasses['Crevasse Number'] == crev_num]
    #             if Df_plot['IsCircular'].all():  ##different colors for circular crevasses and other crevasses
    #                 col = Col_Crevasses_Circ
    #             else:
    #                 continue  ##We don't plot non circular crevasses anymre. Before we had: col = Col_Crevasses_Other
    #             ax8.plot(Df_plot['X'].values / 1000, Df_plot['Y'].values / 1000, color=col, linestyle='-',linewidth=2)
    #             ###Plot cavity contour
    #             ax8.plot(cavity_contour[:, 0] / 1000, cavity_contour[:, 1] / 1000, color=Col_Cavity, linestyle='-',linewidth=2.5)
    #         ###Make a horizontal color bar common to all subplots
    # cbar8 = fig8.colorbar(CS8,
    #                       ax=axes8.ravel().tolist(),  # all subplots
    #                       orientation='horizontal',
    #                       ticks=levs_ticks8,
    #                       fraction=0.04,
    #                       pad=0.08)
    # cbar8.set_label(r'$\sigma_\mathrm{I}$ [kPa]', fontsize=20)
    # cbar8.ax.tick_params(labelsize=18)  # Adjust tick label size
    #
    #
    # ###### START PLOT FIGURE 8 ######
    # fig8bis.suptitle(r'pure Glen transient', fontsize=23, weight='bold')
    # print('starting Figure 8 bis')
    # for i, law in enumerate(['Subhour', 'Subday']):
    #     for j, tsp in enumerate(['a', 'b', 'c', 'd']):
    #         ##get proper ax
    #         ax8 = axes8bis[i, j]
    #         ###Get data for each subplots one by one
    #         df_cav8 = pd.DataFrame()
    #         df_nocav8 = pd.DataFrame()
    #         if i == 0:
    #             if j == 0:  ##Linear elastic
    #                 ax8.set_ylabel(r'Y [km]', fontsize=22)
    #                 df_nocav8 = Data_Simu_El[Data_Simu_El['Case'] == 'NoCavity']
    #                 df_cav8 = Data_Simu_El[Data_Simu_El['Case'] == 'EmptyCavity']
    #                 Title = r'Elastic $\nu = 0.3$'
    #             elif j == 1:  ##Glen after 1s tsp =1s
    #                 df_nocav8 = Data_Simu_VNLTRANS_1s[(Data_Simu_VNLTRANS_1s['Timestep'] == 1) & (Data_Simu_VNLTRANS_1s['Case'] == 'NoCavity')]
    #                 df_cav8 = Data_Simu_VNLTRANS_1s[(Data_Simu_VNLTRANS_1s['Timestep'] == 1) & (Data_Simu_VNLTRANS_1s['Case'] == 'EmptyCavity')]
    #                 Title = r'@1s tsp = 1s'
    #             elif j == 2:  ##Glen after 1mn tsp =1s
    #                 df_nocav8 = Data_Simu_VNLTRANS_1s[(Data_Simu_VNLTRANS_1s['Timestep'] == 60) & (Data_Simu_VNLTRANS_1s['Case'] == 'NoCavity')]
    #                 df_cav8 = Data_Simu_VNLTRANS_1s[(Data_Simu_VNLTRANS_1s['Timestep'] == 60) & (Data_Simu_VNLTRANS_1s['Case'] == 'EmptyCavity')]
    #                 Title = r'@1mn tsp = 1s'
    #             elif j == 3:  ##Glen after 1h tsp=1mn
    #                 df_nocav8 = Data_Simu_VNLTRANS_5mn[(Data_Simu_VNLTRANS_5mn['Timestep'] == 12) & (Data_Simu_VNLTRANS_5mn['Case'] == 'NoCavity')]
    #                 df_cav8 = Data_Simu_VNLTRANS_5mn[(Data_Simu_VNLTRANS_5mn['Timestep'] == 12) & (Data_Simu_VNLTRANS_5mn['Case'] == 'EmptyCavity')]
    #                 Title = r'@1h tsp = 5mn'
    #         elif i == 1:
    #             ax8.set_xlabel(r'X [km]', fontsize=22)
    #             if j == 0:  ##Glen after 6h tsp =5mn
    #                 ax8.set_ylabel(r'Y [km]', fontsize=22)
    #                 df_nocav8 = Data_Simu_VNLTRANS_5mn[(Data_Simu_VNLTRANS_5mn['Timestep'] == 6*60/5) & (Data_Simu_VNLTRANS_5mn['Case'] == 'NoCavity')]
    #                 df_cav8 = Data_Simu_VNLTRANS_5mn[(Data_Simu_VNLTRANS_5mn['Timestep'] == 6*60/5) & (Data_Simu_VNLTRANS_5mn['Case'] == 'EmptyCavity')]
    #                 Title = r'@6h tsp = 5mn'
    #             elif j == 1:  ##Glen after 12h tsp =5mn
    #                 df_nocav8 = Data_Simu_VNLTRANS_5mn[(Data_Simu_VNLTRANS_5mn['Timestep'] == 12*60/5) & (Data_Simu_VNLTRANS_5mn['Case'] == 'NoCavity')]
    #                 df_cav8 = Data_Simu_VNLTRANS_5mn[(Data_Simu_VNLTRANS_5mn['Timestep'] == 12*60/5) & (Data_Simu_VNLTRANS_5mn['Case'] == 'EmptyCavity')]
    #                 Title = r'@12h tsp = 5mn'
    #             elif j == 2:  ##Glen after 24h tsp =5mn
    #                 df_nocav8 = Data_Simu_VNLTRANS_5mn[(Data_Simu_VNLTRANS_5mn['Timestep'] == 24*60/5) & (Data_Simu_VNLTRANS_5mn['Case'] == 'NoCavity')]
    #                 df_cav8 = Data_Simu_VNLTRANS_5mn[(Data_Simu_VNLTRANS_5mn['Timestep'] == 24*60/5) & (Data_Simu_VNLTRANS_5mn['Case'] == 'EmptyCavity')]
    #                 Title = r'@24h tsp = 5mn'
    #             elif j == 3:  ##pure linear viscous incomp
    #                 df_nocav8 = Data_Simu_VNL[Data_Simu_VNL['Case'] == 'NoCavity']
    #                 df_cav8 = Data_Simu_VNL[Data_Simu_VNL['Case'] == 'EmptyCavity']
    #                 Title = r'pure Glen steady'
    #         ###Fill up the plots
    #         ax8.set_title(Title, fontsize=18, weight='bold')
    #         ax8.tick_params(labelsize=18)  # fontsize of the tick labels
    #         ax8.grid(True)
    #         ax8.grid(alpha=0.5)
    #         ###Force x and y axis to same scale
    #         ax8.set_aspect('equal', adjustable='box')
    #         ###Poject SigmaI on interpolation grid
    #         SigmaI_NoCav = Interpolate_field(df_nocav8, 'SigmaI', X, Y)
    #         SigmaI_Cav = Interpolate_field(df_cav8, 'SigmaI', X, Y)
    #         ###shading
    #         clevs8 = clevwide  ## cbar for shading
    #         # colorbar
    #         levs_ticks8 = levs_tickwide
    #         ###Fills up the map with colors for SigmaEq
    #         CS8 = ax8.contourf(X / 1000, Y / 1000, (SigmaI_Cav - SigmaI_NoCav) * 1000, clevs8, cmap=cmap, extend='both')
    #         ax8.plot(xc / 1000, yc / 1000, color='k', linewidth=2)
    #         ###Below we remove colors that are outside of glacier contour
    #         clippath = mpltPath(np.c_[xc / 1000, yc / 1000])
    #         patch = PathPatch(clippath, facecolor='none')
    #         # ax = plt.gca()
    #         ax8.add_patch(patch)
    #         for c in [CS8]:
    #             c.set_clip_path(patch)
    #         ##plot crevasses as continuous line
    #         for crev_num in np.arange(Df_Crevasses['Crevasse Number'].min(),
    #                                   Df_Crevasses['Crevasse Number'].max() + 1):
    #             ###get proper points
    #             Df_plot = Df_Crevasses[Df_Crevasses['Crevasse Number'] == crev_num]
    #             if Df_plot['IsCircular'].all():  ##different colors for circular crevasses and other crevasses
    #                 col = Col_Crevasses_Circ
    #             else:
    #                 continue  ##We don't plot non circular crevasses anymre. Before we had: col = Col_Crevasses_Other
    #             ax8.plot(Df_plot['X'].values / 1000, Df_plot['Y'].values / 1000, color=col, linestyle='-',linewidth=2)
    #             ###Plot cavity contour
    #             ax8.plot(cavity_contour[:, 0] / 1000, cavity_contour[:, 1] / 1000, color=Col_Cavity, linestyle='-',linewidth=2.5)
    #         ###Make a horizontal color bar common to all subplots
    # cbar8bis = fig8bis.colorbar(CS8,
    #                       ax=axes8bis.ravel().tolist(),  # all subplots
    #                       orientation='horizontal',
    #                       ticks=levs_ticks8,
    #                       fraction=0.04,
    #                       pad=0.08)
    # cbar8bis.set_label(r'$\sigma_\mathrm{I}$ [kPa]', fontsize=20)
    # cbar8bis.ax.tick_params(labelsize=18)  # Adjust tick label size



    ###Show map
    plt.show()
