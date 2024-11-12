# -*- coding: utf-8 -*-
"""
Created on Tue Feb  6 12:21:41 2024

@author: dylan

Hodograph rotation using the Taszarek method. This method is done by bring the hodograph down to the 
origin (u0 and v0 = 0 m/s) and rotating the hodograph so the mean wind vector lies at 225 deg.
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
import sharppy.sharptab.utils as utils
import sharppy.sharptab.params as params
import sharppy.sharptab.profile as profile
import sharppy.sharptab.interp as interp
import sharppy.sharptab.winds as winds

warnings.filterwarnings("ignore", category=DeprecationWarning) 
warnings.filterwarnings("ignore", category=FutureWarning)

#plotting the hodograph and the AGL altitude levels
def Hodograph_plot(ax, a, u, v, colour, storm):
    h = Hodograph(ax, component_range=25.)
    h.add_grid(increment=10)
    ax.set_xlabel('Wind Speed (m s\N{SUPERSCRIPT MINUS}\N{SUPERSCRIPT ONE})', size = 12, weight='bold')
    ax.set_ylabel('Wind Speed (m s\N{SUPERSCRIPT MINUS}\N{SUPERSCRIPT ONE})', size = 12, weight='bold')
    
    #finding heights to be labelled on the plot
    hlevs = np.arange(0,12000.0001,10)
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

    h.plot(U_at_hlevs, V_at_hlevs, c=colour, linewidth=2, alpha = 0.6,label=storm)
    ax.legend(loc='upper left')    
    
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

folder_dir = r'C:\Users\dylan\Downloads\NTP Related Documents\2023-2024 Research Projects\NHP Hail Work\SC_unrotated/'

colours=['red', 'blue']

fig = plt.figure(figsize = (8,8))
ax = fig.add_subplot(1, 1, 1)

#iterating over to plot the supercell and non-supercell data using unrotated data (and then rotating)
for j in range(2):
    if j == 0:
        #df = pd.read_csv(folder_dir + 'Hail Profile Averaged.csv')
        df = pd.read_csv(folder_dir + 'SC_RAW_data.csv')
        #storm = 'Supercell'
        storm = 'Hail'
        
        a = df.altitude
        
        u = []
        v = []

        #calculating u and v components (need to give units for metpy)
        for i in range(len(df.ws)):
            u.append((mpcalc.wind_components(df.ws[i] * units('knots').to(units('m/s')),df.wd[i]*units.deg)[0]).magnitude)
            v.append((mpcalc.wind_components(df.ws[i] * units('knots').to(units('m/s')),df.wd[i]*units.deg)[1]).magnitude)
            
        p = df.pressure
        t = df.temp2
        d = df.dpt2

    elif j == 1:
        #df = pd.read_csv(folder_dir + 'NSC_RAW_data.csv')
        #storm = 'Non-Supercell'
        df = pd.read_csv(r'C:\Users\dylan\Downloads\NTP Related Documents\2023-2024 Research Projects\NHP Hail Work\TORNADO_raw_profiles\AB_\Averaged Wind Profile All Months.csv')
        #df = pd.read_csv(folder_dir + 'Region1_AB_RAW_data.csv')
        storm = 'Tornado'
        """
        a = df.altitude
        
        u = []
        v = []

        #calculating u and v components (need to give units for metpy)
        for i in range(len(df.ws)):
            u.append((mpcalc.wind_components(df.ws[i] * units('knots').to(units('m/s')),df.wd[i]*units.deg)[0]).magnitude)
            v.append((mpcalc.wind_components(df.ws[i] * units('knots').to(units('m/s')),df.wd[i]*units.deg)[1]).magnitude)
            
        p = df.pressure
        t = df.temp2
        d = df.dpt2
        """
        u = df.u
        v = df.v
        a = df.alt
        
        p = df.pres
        t = df.temp
        d = df.dpt
        
    #variables needed
    
    
    
    #setting altitude relative to the surface and not sea level
    a = [x - a[0] for x in a]
        
    #setting u0 and v0 to 0 m/s
    u-=u[0]
    v-=v[0] 
    
    ws_new = []
    wd_new = []
    
    #calculating wind speed and direction when u0 and v0 are at 0 m/s
    for i in range(len(u)):
        ws_new.append((mpcalc.wind_speed(u[i]*units('m/s'),v[i]*units('m/s'))).magnitude)
        wd_new.append((mpcalc.wind_direction(u[i]*units('m/s'),v[i]*units('m/s'))).magnitude)
    """
    #creating a profile with sharppy
    prof = profile.create_profile(profile='default', pres=p, hght=a, tmpc=t, \
                                        dwpc=d, wspd=ws_new, wdir=wd_new, missing=-9999, strictQC=True)
    sfc = prof.pres[prof.sfc]
    p6km = interp.pres(prof, interp.to_msl(prof, 6000.))

    # Calculate the 0-6 km pressure-weighted mean wind
    mean_6km = winds.mean_wind(prof, pbot=sfc, ptop=p6km)
    mean_6kmu = mean_6km[0]
    mean_6kmv = mean_6km[1]
    """
    hlevs = np.arange(0,6000.0001,10)
    U_interpf = interp1d(a, u, kind='linear', bounds_error=False)
    V_interpf = interp1d(a, v, kind='linear', bounds_error=False)
    U_at_hlevs = U_interpf(hlevs)
    V_at_hlevs = V_interpf(hlevs)
    u0=(mean(U_at_hlevs))
    v0=(mean(V_at_hlevs))
    
    #finding the mean wind direction
    wd_mean = mpcalc.wind_direction(u0*units('m/s'),v0*units('m/s')).magnitude
    
    #wind speed and wind direction of the 0-6 km mean wind
    #ws_mean = mpcalc.wind_speed(mean_6kmu*units['m/s'], mean_6kmv*units['m/s'])
    #wd_mean = (mpcalc.wind_direction(mean_6kmu*units['m/s'], mean_6kmv*units['m/s'])).magnitude
    
    #rotating hodograph so the mean wind is at 225 deg
    rot = 225-wd_mean
    wd_new_rot = wd_new+rot
    
    #calculating new u and v components to be used for the hodograph
    u_new = []
    v_new = []
    
    #calculating the new u and v components of the wind
    for i in range(len(ws_new)):
        u_new.append((mpcalc.wind_components(ws_new[i] * units('m/s'),wd_new_rot[i]*units.deg)[0]).magnitude)
        v_new.append((mpcalc.wind_components(ws_new[i] * units('m/s'),wd_new_rot[i]*units.deg)[1]).magnitude)
    
    #plotting the hodograph
    Hodograph_plot(ax,a,u_new,v_new,colours[j], storm)
    #Hodograph_plot(ax,a,u,v,colours[j], storm)
    

#folder_dir = r'C:\Users\dylan\Downloads\NHP Hail Work\SC_rotated/'
"""
#plotting already rotated data
colours = ['cyan', 'green']
for j in range(2):
    if j == 0:
        df = pd.read_csv(folder_dir + 'SC_RAW_data.csv')
        #storm = 'Supercell'
        storm = 'Hail'
    elif j == 1:
        #df = pd.read_csv(folder_dir + 'NSC_RAW_data.csv')
        #storm = 'Non-Supercell'
        df = pd.read_csv(folder_dir + 'Region1_AB_RAW_data.csv')
        storm = 'Tornado'
        
    #variables needed
    a = df.altitude
    p = df.pressure
    t = df.temp2
    d = df.dpt2
    ws = df.ws
    wd = df.wd
    
    #setting altitude to AGL
    a = [x - a[0] for x in a]

    u = []
    v = []

    #finding u and v components of the wind (metpy needs units for calculation)
    for i in range(len(df.ws)):
        u.append((mpcalc.wind_components(df.ws[i] * units('m/s'),df.wd[i]*units.deg)[0]).magnitude)
        v.append((mpcalc.wind_components(df.ws[i] * units('m/s'),df.wd[i]*units.deg)[1]).magnitude)
    
    #plotting hodograph
    Hodograph_plot(ax,a,u,v,colours[j], storm)
    
"""
title = 'Average Hail - Tornado Supercell Rotated - Taszarek Method'
#title = 'Average Hail - Tornado AB Supercell Unrotated'
ax.set_title(f'{title}')
#plt.savefig(r'C:\Users\dylan\Downloads\NHP Hail Work\Plots\Correct Plots/'+f'{title} Hodograph.png',dpi=300, bbox_inches='tight')
#plt.savefig(r'C:\Users\dylan\Downloads\NHP Hail Work\Plots\Correct Plots/'+f'{title} Hodograph.png',dpi=300, bbox_inches='tight')