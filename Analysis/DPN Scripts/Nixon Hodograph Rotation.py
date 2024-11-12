# -*- coding: utf-8 -*-
"""
Created on Tue Feb  6 12:21:41 2024

@author: dylan

Scipt looking at rotating the hodograph using the Nixon-Allen (2022) approach. This is done rotating
the hodograph so the 0-3km shear vector lies along the positive x-axis. The hodograph is then brought
to the orgin - so u0 and v0 = 0m/s.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from statistics import mean
from metpy.plots import Hodograph
import metpy.calc as mpcalc
from metpy.units import units
from scipy.interpolate import interp1d
import warnings
from matplotlib.lines import Line2D

warnings.filterwarnings("ignore", category=DeprecationWarning) 
warnings.filterwarnings("ignore", category=FutureWarning)

#plotting the hodograph and the AGL altitude levels
def Hodograph_plot(ax, a, u, v, colour, title, storm):
    h = Hodograph(ax, component_range=30.)
    h.add_grid(increment=10)
    ax.set_xlabel('Wind Speed (m s\N{SUPERSCRIPT MINUS}\N{SUPERSCRIPT ONE})', size = 12, weight='bold')
    ax.set_ylabel('Wind Speed (m s\N{SUPERSCRIPT MINUS}\N{SUPERSCRIPT ONE})', size = 12, weight='bold')
    
    #finding heights to be labelled on the plot
    hlevs = np.arange(0,12200.0001,100)
    lev0km_ind = np.where(hlevs==0)[0]
    lev1km_ind = np.where(hlevs==1000)[0]
    lev3km_ind = np.where(hlevs==3000)[0]
    lev6km_ind = np.where(hlevs==6000)[0]
    lev9km_ind = np.where(hlevs==9000)[0]
    lev12km_ind = np.where(hlevs==12000)[0]

    #interpolating heights
    U_interpf = interp1d(a, u, kind='linear', bounds_error=False)
    V_interpf = interp1d(a, v, kind='linear', bounds_error=False)
    U_at_hlevs = U_interpf(hlevs)
    V_at_hlevs = V_interpf(hlevs)

    us = []
    vs = []
    ws = []
    wd = []
    wd2 = []

    
    #calculating wind speeds and directions at each level so they can be used to rotate the hodograph
    for i in range(len(U_at_hlevs)):
        ws.append(mpcalc.wind_speed(U_at_hlevs[i]*units('m/s'),V_at_hlevs[i]*units('m/s')))
        wd.append(mpcalc.wind_direction(U_at_hlevs[i]*units('m/s'),V_at_hlevs[i]*units('m/s')))
        
    #rotating hodogaph so the 0-3km shear vector lies along the postive x-axis, before calculating the
    #u and v components (so hodograph can be brought down to the origin)
    for i in range(len(ws)):
        wd2.append(wd[i].magnitude + (270-(wd[int(lev3km_ind)].magnitude)))
        #wd2.append(wd[i].magnitude)
        us.append((mpcalc.wind_components(ws[i],wd2[i]*units.deg)[0]).magnitude)
        vs.append((mpcalc.wind_components(ws[i],wd2[i]*units.deg)[1]).magnitude)
    
    #bringing hodograph down to the origin
    us-=us[0]
    vs-=vs[0] 
    
    #plotting of hodograph
    h.plot(us, vs, c=colour, linewidth=2, alpha = 0.6, label=storm)
    ax.legend(loc='upper left')    
    
    #adding text to plot
    ax.text(us[int(lev0km_ind)],vs[int(lev0km_ind)],str(int(hlevs[int(lev0km_ind)]/1000)),color='k',
             fontsize=9,horizontalalignment='center',verticalalignment='center',
             bbox=dict(boxstyle='circle',facecolor='gray',edgecolor='none',alpha=0.3,pad=0.1))
    ax.text(us[int(lev1km_ind)],vs[int(lev1km_ind)],str(int(hlevs[int(lev1km_ind)]/1000)),color='k',
             fontsize=9,horizontalalignment='center',verticalalignment='center',
             bbox=dict(boxstyle='circle',facecolor='gray',edgecolor='none',alpha=0.3,pad=0.1))
    ax.text(us[int(lev3km_ind)],vs[int(lev3km_ind)],str(int(hlevs[int(lev3km_ind)]/1000)),color='k',
             fontsize=9,horizontalalignment='center',verticalalignment='center',
             bbox=dict(boxstyle='circle',facecolor='gray',edgecolor='none',alpha=0.3,pad=0.1))
    ax.text(us[int(lev6km_ind)],vs[int(lev6km_ind)],str(int(hlevs[int(lev6km_ind)]/1000)),color='k',
             fontsize=9,horizontalalignment='center',verticalalignment='center',
             bbox=dict(boxstyle='circle',facecolor='gray',edgecolor='none',alpha=0.3,pad=0.1))
    ax.text(us[int(lev9km_ind)],vs[int(lev9km_ind)],str(int(hlevs[int(lev9km_ind)]/1000)),color='k',
             fontsize=9,horizontalalignment='center',verticalalignment='center',
             bbox=dict(boxstyle='circle',facecolor='gray',edgecolor='none',alpha=0.3,pad=0.1))
    ax.text(us[int(lev12km_ind)],vs[int(lev12km_ind)],str(int(hlevs[int(lev12km_ind)]/1000)),color='k',
             fontsize=9,horizontalalignment='center',verticalalignment='center',
             bbox=dict(boxstyle='circle',facecolor='gray',edgecolor='none',alpha=0.3,pad=0.1))
    
plt.close('all')

folder_dir = r'C:\Users\dylan\Downloads\NTP Related Documents\2023-2024 Research Projects\NHP Hail Work\SC_unrotated/'

colours=['red', 'blue']

fig = plt.figure(figsize = (8,8))
ax = fig.add_subplot(1, 1, 1)

title = 'Rotated Hodographs'

#iterating over supercell (supercell hail) and non-supercell (supercell tornado) events
for j in range(2):
    if j == 0:
        #df = pd.read_csv(folder_dir + 'Hail Profile Averaged.csv')
        df = pd.read_csv(folder_dir + 'SC_RAW_data.csv')
        storm = 'Supercell'
        #storm = 'Hail'
        
        a = df.altitude
        u = []
        v = []

        #calculating u and v components (need to give units for metpy)
        for i in range(len(df.ws)):
            u.append((mpcalc.wind_components(df.ws[i] * units('knots').to(units('m/s')),df.wd[i]*units.deg)[0]).magnitude)
            v.append((mpcalc.wind_components(df.ws[i] * units('knots').to(units('m/s')),df.wd[i]*units.deg)[1]).magnitude)
        
    elif j == 1:
        df = pd.read_csv(folder_dir + 'NSC_RAW_data.csv')
        storm = 'Non-Supercell'
        #df = pd.read_csv(r'C:\Users\dylan\Downloads\NHP Hail Work\TORNADO_raw_profiles\AB_\Averaged Wind Profile All Months.csv')
        #df = pd.read_csv(folder_dir + 'Region1_AB_RAW_data.csv')
        #storm = 'Tornado'
        """
        a = df.alt
        u = df.u
        v = df.v
        """
        a = df.altitude
        u = []
        v = []

        #calculating u and v components (need to give units for metpy)
        for i in range(len(df.ws)):
            u.append((mpcalc.wind_components(df.ws[i] * units('knots').to(units('m/s')),df.wd[i]*units.deg)[0]).magnitude)
            v.append((mpcalc.wind_components(df.ws[i] * units('knots').to(units('m/s')),df.wd[i]*units.deg)[1]).magnitude)
    """
    u = []
    v = []
    for i in range(len(df.ws)):
        u.append((mpcalc.wind_components(df.ws[i] * units('knots').to(units('m/s')),df.wd[i]*units.deg)[0]).magnitude)
        v.append((mpcalc.wind_components(df.ws[i] * units('knots').to(units('m/s')),df.wd[i]*units.deg)[1]).magnitude)
    """
    
    a = [x - a[0] for x in a]
    
    Hodograph_plot(ax,a,u,v,colours[j],title, storm)


title = 'Average Hail SC -NSC Rotated - Nixon & Allen Method'
#title = 'Average Supercell - Non-Supercell Nixon & Allen Rotated'
ax.set_title(f'{title}')
#plt.savefig(r'C:\Users\dylan\Downloads\NHP Hail Work\Plots\Correct Plots/' + f'{title} Hodograph-Updated.png',dpi=300, bbox_inches='tight')