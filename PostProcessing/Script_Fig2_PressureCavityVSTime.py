################################################################################
# import #####
################################################################################
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cmx

from pathlib import Path

import numpy as np
from matplotlib.ticker import FormatStrFormatter
import matplotlib.gridspec as gridspec

import pandas as pd  ###To treat the csv files
from datetime import datetime, date, timedelta ###To handle dates
from matplotlib.dates import MonthLocator, DayLocator, HourLocator, DateFormatter, drange
from matplotlib import ticker
from scipy import interpolate

# formatter = ticker.ScalarFormatter(useMathText=True)
# formatter.set_scientific(True)
# formatter.set_powerlimits((-3, 3))

###Defined Colors optimized for color-blind people:
RedStar = [162 / 255, 20 / 255, 47 / 255]
RedBackGround = [233 / 255, 201 / 255, 207 / 255]
LightBlueGreyBackGround = [194 / 255, 209 / 255, 219 / 255]
BlueGreyBackGround = [147/ 255, 174 / 255, 191 / 255]
LightBrownBackGround = [223/ 255, 211 / 255, 205 / 255]
VeryLightBrown = [255 / 255, 199 / 255, 127 / 255]
LightBrown = [229 / 255, 143 / 255, 91 / 255]
Brown = [140 / 255, 88 / 255, 56 / 255]
DarkBrown = [64 / 255, 40 / 255, 25 / 255]
VeryDarkBrown = [13 / 255, 8 / 255, 5 / 255]

####//////////////////////////////////////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\####
####///~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\\\####
####///~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ MAIN PART OF THE CODE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\\\####
####///~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\\\####
####//////////////////////////////////////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\####

################################################################################
# prepare plots #####
################################################################################
# global figure pimping
# plt.rc('text', usetex=True)
plt.rc('xtick', labelsize=22)
plt.rc('ytick', labelsize=22)
plt.rc('axes', labelsize=22)
plt.rc('legend', fontsize=26)

###Fig 1 is Water level as a funciton of time
fig1 = plt.figure(1, figsize=(40, 20))
plt.ylabel(r'Water Level (m)', fontsize=28)
# plt.xlabel(r'Day of Simu', fontsize=34)
# plt.xlim([0, 65])
plt.tick_params(labelsize=24)  # fontsize of the tick labels
plt.grid(True)

###############################################################################
# Make the plots #####
################################################################################
if __name__ == "__main__":
    RefStartDate = datetime(2011, 7, 1) ##Every days are numbered relative to 1st July 2011 (day 1)
    SimuStartDay = -315 ###Number of the day corresponding to the start of the simu pumping 2010
    Pathroot = Path('/home/brondexj/BETTIK/TeteRousse/MyTeterousse_GeoGag/Data/.')
    Filename = 'Evol_Niveau2010-2013-moyen.dat'
    Col_Names = ['Day', 'WaterLevel']
    Df_WaterLevel = pd.read_csv(Pathroot.joinpath(Filename), names=Col_Names, delim_whitespace=True, skiprows=1)
    ### Create a column converting day number into date
    Df_WaterLevel['Date'] = RefStartDate + pd.to_timedelta(Df_WaterLevel['Day'], unit='D')
    ### Simu Pumping 2010 start day one and give one output every day:
    DayOfSimuPumping2010 = np.arange(SimuStartDay, SimuStartDay+61, 1)
    ### Values of Water Level on Simu days are interpolated from water level measured every now and then (Data file)
    WaterLevel_Interp=np.interp(DayOfSimuPumping2010,Df_WaterLevel['Day'],Df_WaterLevel['WaterLevel']) ##Linear Interpolation as in Elmer
    Colormap_for_Simudays = 'copper_r'

    ###ColorMap for correspondance between Water Pressure and stress
    fig1 = plt.figure(1)
    plt.tick_params(labelsize=22)  # fontsize of the tick labels
    # Set xlim from first day of simu (20 August 2010) to last day of simu (05 September 2013)
    plt.xlim(pd.to_datetime('2010-08-20'), pd.to_datetime('2013-09-05'))
    plt.gca().xaxis.set_major_locator(MonthLocator(interval=2))
    plt.gca().xaxis.set_minor_locator(MonthLocator(interval=1))
    plt.gca().xaxis.set_major_formatter(DateFormatter('%Y-%m'))
    fig1.autofmt_xdate()
    ###Shade periods corresponding to pumping (26/08/2010 to 15/10/10; 28/09/11 to14/10/11; 23/09/12 to 09/10/12)
    plt.axvspan(date(int(2010), 8, 26), date(int(2010), 10, 15), alpha=0.3, color='grey')
    plt.axvspan(date(int(2011), 9, 28), date(int(2011), 10, 14), alpha=0.3, color='grey')
    plt.axvspan(date(int(2012), 9, 23), date(int(2012), 10, 9), alpha=0.3, color='grey')
    ###Plot the pressure over the full period
    plt.plot_date( Df_WaterLevel['Date'].values, Df_WaterLevel['WaterLevel'].values, color='k', linestyle='-', linewidth=1.7, marker='None')

    ### Here plot color dots for the days at which we have an equivalent stress profile

    cm = plt.get_cmap(Colormap_for_Simudays)
    cNorm = colors.Normalize(vmin=np.min(DayOfSimuPumping2010), vmax=np.max(DayOfSimuPumping2010))
    scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=cm)
    # plt.scatter(1, WaterLevel_Interp[0], facecolor=scalarMap.to_rgba(1), edgecolors = 'k', s=600, zorder=3)
    for j, Dday in enumerate(DayOfSimuPumping2010):
        ### for pumping step profile every 5 days
        if not Dday % 5 == 0:
            continue
        ### Convert day in date
        CorrespDate = RefStartDate+timedelta(days=int(Dday))
        ###Scatter Plot of Water Level for considered day
        if not Dday == np.max(DayOfSimuPumping2010):
            plt.scatter(CorrespDate, WaterLevel_Interp[j], facecolor=scalarMap.to_rgba(Dday), edgecolors = 'k', s=100, zorder=3)
        else:
            plt.scatter(CorrespDate, WaterLevel_Interp[j], facecolor=scalarMap.to_rgba(Dday), marker='d', edgecolors = 'k', s=150, zorder=3)
    plt.show()