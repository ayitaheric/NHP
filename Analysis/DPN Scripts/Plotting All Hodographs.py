# -*- coding: utf-8 -*-
"""
Created on Tue Feb  6 12:21:41 2024

@author: dylan

Plots each hodograph within the data and plots the Bunkers RM and LM
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from statistics import mean
from metpy.plots import Hodograph
import metpy.calc as mpcalc
from metpy.units import units, pandas_dataframe_to_unit_arrays
from scipy.interpolate import interp1d
import warnings
from matplotlib.lines import Line2D

warnings.filterwarnings("ignore", category=DeprecationWarning) 
warnings.filterwarnings("ignore", category=FutureWarning)

#hodograph plot
def Hodograph_plot(ax, a, u, v, colour):#, title):
    h = Hodograph(ax, component_range=55.)
    h.add_grid(increment=10)
    ax.set_xlabel('Wind Speed (m s\N{SUPERSCRIPT MINUS}\N{SUPERSCRIPT ONE})', size = 12, weight='bold')
    ax.set_ylabel('Wind Speed (m s\N{SUPERSCRIPT MINUS}\N{SUPERSCRIPT ONE})', size = 12, weight='bold')
    
    hlevs = np.arange(0,12200.0001,200)
    lev0km_ind = np.where(hlevs==0)[0]
    lev1km_ind = np.where(hlevs==1000)[0]
    lev3km_ind = np.where(hlevs==3000)[0]
    lev6km_ind = np.where(hlevs==6000)[0]
    lev9km_ind = np.where(hlevs==9000)[0]
    lev12km_ind = np.where(hlevs==12000)[0]

    U_interpf = interp1d(a, u, kind='linear', bounds_error=False)
    V_interpf = interp1d(a, v, kind='linear', bounds_error=False)
    U_at_hlevs = U_interpf(hlevs)
    V_at_hlevs = V_interpf(hlevs)

    h.plot(U_at_hlevs, V_at_hlevs, c=colour, linewidth=2, alpha = 0.6)
    
    ax.text(U_at_hlevs[lev0km_ind],V_at_hlevs[lev0km_ind],str(int(hlevs[lev0km_ind]/1000)),color='k',
             fontsize=9,horizontalalignment='center',verticalalignment='center',
             bbox=dict(boxstyle='circle',facecolor='gray',edgecolor='none',alpha=0.3,pad=0.1))
    ax.text(U_at_hlevs[lev1km_ind],V_at_hlevs[lev1km_ind],str(int(hlevs[lev1km_ind]/1000)),color='k',
             fontsize=9,horizontalalignment='center',verticalalignment='center',
             bbox=dict(boxstyle='circle',facecolor='gray',edgecolor='none',alpha=0.3,pad=0.1))
    ax.text(U_at_hlevs[lev3km_ind],V_at_hlevs[lev3km_ind],str(int(hlevs[lev3km_ind]/1000)),color='k',
             fontsize=9,horizontalalignment='center',verticalalignment='center',
             bbox=dict(boxstyle='circle',facecolor='gray',edgecolor='none',alpha=0.3,pad=0.1))
    ax.text(U_at_hlevs[lev6km_ind],V_at_hlevs[lev6km_ind],str(int(hlevs[lev6km_ind]/1000)),color='k',
             fontsize=9,horizontalalignment='center',verticalalignment='center',
             bbox=dict(boxstyle='circle',facecolor='gray',edgecolor='none',alpha=0.3,pad=0.1))
    ax.text(U_at_hlevs[lev9km_ind],V_at_hlevs[lev9km_ind],str(int(hlevs[lev9km_ind]/1000)),color='k',
             fontsize=9,horizontalalignment='center',verticalalignment='center',
             bbox=dict(boxstyle='circle',facecolor='gray',edgecolor='none',alpha=0.3,pad=0.1))
    ax.text(U_at_hlevs[lev12km_ind],V_at_hlevs[lev12km_ind],str(int(hlevs[lev12km_ind]/1000)),color='k',
             fontsize=9,horizontalalignment='center',verticalalignment='center',
             bbox=dict(boxstyle='circle',facecolor='gray',edgecolor='none',alpha=0.3,pad=0.1))
    
plt.close('all')

folder_dir = r'C:\Users\dylan\Downloads\NTP Related Documents\2023-2024 Research Projects\NHP Hail Work\Hail_data_MATEUSZ\Raw_profiles/'

df_reports = pd.read_csv(r'C:\Users\dylan\Downloads\NTP Related Documents\2023-2024 Research Projects\NHP Hail Work\Hail_data_MATEUSZ\Chosen Profiles.csv')
df_reports['datetime'] = pd.to_datetime(df_reports['ERA5_DATE'])

hlevs = np.arange(0,12200.0001,200)
lev0km_ind = np.where(hlevs==0)[0]
lev1km_ind = np.where(hlevs==1000)[0]
lev3km_ind = np.where(hlevs==3000)[0]
lev6km_ind = np.where(hlevs==6000)[0]
lev9km_ind = np.where(hlevs==9000)[0]
lev12km_ind = np.where(hlevs==12000)[0]

#iterating over each profile to plot the hodographs
for i in range(len(df_reports['Event Number'])):
    fig = plt.figure(figsize = (8,8))
    ax = fig.add_subplot(1, 1, 1)
    
    time = df_reports.datetime[i].strftime('_%Y-%m-%d_')
    
    if int(df_reports['Event Number'][i]) <10:
        num = f'_000{df_reports["Event Number"][i]}'
    else:
        num = f'_00{df_reports["Event Number"][i]}'
        
    hour = f'{df_reports.datetime[i].strftime("%H")}_'
        
    if df_reports['Storm Type'][i] == 'Y':
        colours = 'red'
        storm = 'Supercell'
    else:
        colours = 'blue'
        storm = 'Non-Supercell'
    
    df = pd.read_csv(folder_dir + f'HAIL{num}{time}{hour}UTC.csv')
        
    a = df.altitude
    p = df.pressure
    u = df.u
    v = df.v
    a = [x - a[0] for x in a] #seeting altitude from ASL to AGL
    
    Hodograph_plot(ax,a,u,v,colours[0])
    
    unit_types = {'pressure':'hPa', 'temp':'degC', 'dpt':'degC', 'wd':'deg', 'ws':'m/s',\
                  'u':'m/s', 'v':'m/s'}
    united_array=pandas_dataframe_to_unit_arrays(df,column_units=unit_types)

    a = united_array['altitude']
    p = united_array['pressure']
    u = united_array['u']
    v = united_array['v']
    brm = (mpcalc.bunkers_storm_motion(p, u, v, height = a*units.m)[0]).magnitude
    blm = (mpcalc.bunkers_storm_motion(p, u, v, height = a*units.m)[1]).magnitude
    
    #plotting BRM and BLM
    ax.plot(brm[0], brm[1], marker='o', markersize=8, markerfacecolor = 'red', markeredgecolor = 'black')
    ax.plot(blm[0], blm[1], marker='o', markersize=8, markerfacecolor = 'blue', markeredgecolor = 'black')
    
    title = f'{df_reports.ERA5_DATE[i]} - {storm}'
    fig.suptitle(f'{title}')

    #plt.savefig(r'C:\Users\dylan\Downloads\NHP Hail Work\Plots\All Hodographs/' + f'{df_reports.datetime[i].strftime("%Y_%m_%d_%H")}.png',dpi=300, bbox_inches='tight')