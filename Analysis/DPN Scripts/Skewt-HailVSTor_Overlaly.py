# -*- coding: utf-8 -*-
"""
Created on Tue Feb  6 12:21:41 2024

@author: dylan

Script where you can overlay two (or more if desired) skew-t plots. ML LCL, LFC, EL, CAPE and CIN 
(shaded) are included.
"""

import pandas as pd
from datetime import timedelta
import numpy as np
import matplotlib.pyplot as plt
from statistics import mean
from metpy.plots import SkewT, Hodograph
import metpy.calc as mpcalc
from metpy.units import units
import matplotlib.colors as mcolors
import matplotlib.patheffects as mpatheffects
import metpy.interpolate as mpinterpolate
from scipy.interpolate import interp1d
import warnings
import sharppy.sharptab.profile as profile
import sharppy.sharptab.interp as interp
import sharppy.sharptab.winds as winds
import sharppy.sharptab.utils as utils
import sharppy.sharptab.params as params
import sharppy.sharptab.thermo as thermo
import wrf

warnings.filterwarnings("ignore", category=DeprecationWarning) 
warnings.filterwarnings("ignore", category=FutureWarning)

plt.close('all')

folder_dir = r'C:\Users\dylan\Downloads\NTP Related Documents\2023-2024 Research Projects\NHP Hail Work\SC_rotated/'
#arr = np.empty((5,), dtype=np.ndarray)

fig = plt.figure(figsize=(10,10))
skew = [np.nan]*2
skew = SkewT(fig, subplot = (1,1,1))

temps = ['red', 'darkorange']
dews = ['green', 'blue']
tvs = ['fuchsia', 'purple']
capes = ['lightcoral', 'gold']
cins = ['cyan', 'dodgerblue']

#iterating over the data for the two skew-ts
for h in range(2):
    if h ==0:
        df = pd.read_csv(folder_dir + 'SC_RAW_data.csv')
        storm = 'Hail'
        #storm = 'SC'
    elif h ==1:
        df = pd.read_csv(folder_dir + 'Region1_AB_RAW_data.csv')
        storm = 'Tornado'
        #Region1_AB_RAW_data
        #df = pd.read_csv(folder_dir + 'NSC_RAW_data.csv')
        #storm = 'NSC'
        
    a = df.altitude
    p = df.pressure
    t = df.temp2
    d = df.dpt2
    
    #calculating u and v components for calculation
    u = []
    v = []
    data_tv = []
    for i in range(len(df.ws)):
        u.append((mpcalc.wind_components(df.ws[i] * units('knots').to(units('m/s')),df.wd[i]*units.deg)[0]).magnitude)
        v.append((mpcalc.wind_components(df.ws[i] * units('knots').to(units('m/s')),df.wd[i]*units.deg)[1]).magnitude)
        data_tv.append(thermo.virtemp(p[i], t[i], d[i]))
        
    data_tv = np.array(data_tv)
    
    #setting altitude to AGL instead of ASL
    a = [x - a[0] for x in a]

    wdir = mpcalc.wind_direction(u*units('m/s'),v*units('m/s'))
    wspd = mpcalc.wind_speed(u*units('m/s'),v*units('m/s'))
    
    #creating profile
    data_profset = profile.create_profile(profile='default', pres=p, hght=a,
                                          tmpc=t, dwpc=d, wdir=wdir, wspd=wspd)

    sfcpres = data_profset.pres[data_profset.sfc]
    # sfcpres = profset.pres[profset.sfc]

    # p3km = interp.pres(prof, interp.to_msl(prof, 3000.))
    
    #finding values related to the parcel
    parcel_sfcp = data_profset.pres[0]
    parcel_topp = interp.pres(data_profset, interp.to_msl(data_profset, 500.))
    parcel_depth = parcel_sfcp - parcel_topp
    
    # Define various parcels to lift
    data_sfcpcl = params.parcelx(data_profset, flag=1)  # Surface-based parcel
    data_fcstpcl = params.parcelx(data_profset, flag=2) # Forecast parcel (not sure what this means, though)
    data_mupcl = params.parcelx(data_profset, flag=3)   # Most-unstable parcel
    # alternative workaround
    mlpcl_vals = params.DefineParcel(data_profset, flag=4, pres=parcel_depth)
    # mlpcl_vals = params.DefineParcel(data_profset, flag=4, pres=100.)
    data_mlpcl = params.parcelx(data_profset, lplvals=mlpcl_vals)
    
    #plotting data
    skew.plot(p,t,c=temps[h],linewidth=1,ls='-',label=f'Temperature - {storm}')
    skew.plot(p,data_tv, c=temps[h], linewidth=1,ls=':',label=f'Virtual Temperature - {storm}')
    skew.plot(p,d,c=dews[h],linewidth=1,ls='-',label=f'Dewpoint - {storm}')
    # skews.plot(data_sfcpcl.ptrace,data_sfcpcl.ttrace,color=[0.5,0,0.5],linewidth=1, linestyle='--')
    skew.plot(data_mlpcl.ptrace,data_mlpcl.ttrace,color=tvs[h],linewidth=1, linestyle='-',label=f'ML Parcel Trace - {storm}')
    #skew.plot_barbs(p[0:30:5],u[0:30:5],v[0:30:5],color='k',sizes={'spacing':0.2},length=5.,xloc=1.09)
    #skew.plot_barbs(p[30::3],u[30::3],v[30::3],color='k',sizes={'spacing':0.2},length=5.,xloc=1.09)

    parc_tvshade = wrf.interp1d(np.round(data_mlpcl.ttrace.data,4), np.round(data_mlpcl.ptrace.data,4), np.round(p,5))
    # parc_tdshade = wrf.interp1d(np.round(data_mlpcl.ptrace.data,5), np.round(data_mlpcl.ptrace.data,5), np.round(data_pres,5))
    parc_pshade = wrf.interp1d(np.round(data_mlpcl.ptrace.data,4), np.round(data_mlpcl.ptrace.data,4), np.round(p,5))
    
    parc_tvshade[0] = data_mlpcl.ttrace.data[0]
    parc_pshade[0] = data_mlpcl.ptrace.data[0]
    
    # env_tvshade = wrf.interp1d(data_tv, data_pres, np.round(data_mlpcl.ptrace.data,6))
    # env_pshade = wrf.interp1d(data_pres, data_pres, np.round(data_mlpcl.ptrace.data,6))
    # env_tdshade = wrf.interp1d(data_dew, data_pres, np.round(data_mlpcl.ptrace.data,6))

    # envdew_shade = wrf.interp1d(data_dew, data_pres, data_mlpcl.ptrace)
    
    lcl_pres = data_mlpcl.lclpres
    lcl_temp = wrf.interp1d(np.array(data_tv), np.array(p), np.asarray([lcl_pres]))
    # print('Region:',region[fin], "ML LCL pressure:", lcl_pres)
    lfc_pres = data_mlpcl.lfcpres
    lfc_temp = wrf.interp1d(np.array(data_tv), np.array(p), np.asarray([lfc_pres]))
    el_pres = data_mlpcl.elpres
    el_temp = wrf.interp1d(np.array(data_tv), np.array(p), np.asarray([el_pres]))
    parc_eltemp = wrf.interp1d(parc_tvshade, parc_pshade, np.asarray([el_pres]))
    
    # skews[fin].plot(data_mlpcl.pres, data_mlpcl.tmpc, 'ro', markersize=2, markerfacecolor='red')
    
    cap = data_mlpcl.cap
    ccd = data_mlpcl.elhght-data_mlpcl.lfchght
    
    # sys.exit()
    
    if h == 0:
        align = 'left'
        value = 1.5
    else:
        align = 'right'
        value = -1.5
    
    #plotting ML LCL, LFC, EL
    skew.plot([lcl_pres, lcl_pres], [lcl_temp-1, lcl_temp+1], 'k-', lw=1)
    skew.ax.text(lcl_temp.data+value, lcl_pres, f'ML LCL - {storm}',fontsize=8,verticalalignment='center',horizontalalignment=f'{align}',bbox=dict(facecolor='gray',edgecolor='none',alpha=0.5,pad=0.1))
    skew.plot([lfc_pres, lfc_pres], [lfc_temp-1, lfc_temp+1], 'k-', lw=1)
    skew.ax.text(lfc_temp.data+value, lfc_pres, f'ML LFC - {storm}',fontsize=8,verticalalignment='center',horizontalalignment=f'{align}',bbox=dict(facecolor='gray',edgecolor='none',alpha=0.5,pad=0.1))
    skew.plot([el_pres, el_pres], [el_temp-1, el_temp+1], 'k-', lw=1)
    skew.ax.text(el_temp.data+value, el_pres, f'ML EL - {storm}',fontsize=8,verticalalignment='center',horizontalalignment=f'{align}',bbox=dict(facecolor='gray',edgecolor='none',alpha=0.5,pad=0.1))
    
    # skews[fin].shade_cin(data_pshade., data_tshade, data_mlpcl.ttrace.data, dewpoint=data_tdshade)
    # skews[fin].shade_cape(data_pshade.data, data_tshade, data_mlpcl.ttrace.data, dewpoint=data_tdshade)

    # skews[fin].shade_cin(data_mlpcl.ptrace, envtempv_shade, data_mlpcl.ttrace.data)
    # skews[fin].shade_cape(data_mlpcl.ptrace, envtempv_shade, data_mlpcl.ttrace.data)

    # find the level where environment Tv exceeds the parcel TV again.
    EL_lvlind = np.min(np.where(p < el_pres))-1
    
    data_pres_cinshade = np.append(p[0:EL_lvlind], el_pres)
    data_temp_cinshade = np.append(data_tv[0:EL_lvlind], el_temp)  # really the environmental virtual temperature
    parc_tv_cinshade = np.append(parc_tvshade[0:EL_lvlind], parc_eltemp)

    # skews[fin].shade_cin(data_pres[0:EL_lvlind], data_tv[0:EL_lvlind], parc_tvshade[0:EL_lvlind])
    skew.shade_cin(data_pres_cinshade, data_temp_cinshade, parc_tv_cinshade,color=cins[h],alpha=0.3,label=f'CIN - {storm}')
    skew.shade_cape(p, data_tv, parc_tvshade, alpha = 0.2,color=capes[h],label=f'CAPE - {storm}')
    """
    def pressure_interval(p,u,v,upper=100,lower=1000,spacing=50):
    
        intervals = list(range(upper,lower,spacing))
    
        ix = []
        for center in intervals:
            index = (np.abs(np.array(p)-center)).argmin()
            if index not in ix:
                ix.append(index)
                
        ps = []
        us = []
        vs = []
        for i in ix:
            ps.append(p[i])
            us.append(u[i])
            vs.append(v[i])
    
        return ps,us,vs
    
    p_,u_,v_ = pressure_interval(p,u,v)
    skew[h].plot_barbs(p_,u_,v_)
    """
    
    skew.ax.set_xlim(-40,50)
    skew.ax.set_ylim(1000,100)
    skew.ax.set_xticks(np.arange(-40,50.01,10))
    skew.ax.set_yticks(np.arange(1000,99.99,-100))
    skew.ax.set_xticklabels(['-40','-30','-20','-10','0','10','20','30','40','50'],fontsize=11)
    skew.ax.set_yticklabels(['1000','900','800','700','600','500','400','300','200','100'],fontsize=11)
    skew.ax.set_xlabel('Temperature ($\mathregular{^o}$C)',fontsize=12)
    skew.ax.set_ylabel('Pressure (hPa)', fontsize=12)
    #skew.ax.axvline(0,color='c',linestyle='--',linewidth=1.5)
    skew.plot_dry_adiabats(t0=(np.arange(-50,260,10))*units('degC'),color='gray',linewidth=0.5,linestyle='--', alpha = 0.3)
    skew.plot_moist_adiabats(color='gray',linewidth=0.5,linestyle='--', alpha = 0.3)
    skew.plot_mixing_lines(color='gray',linewidth=0.5,linestyle=':')
    skew.ax.axvline(-30, color='darkblue', linestyle='--', alpha=0.7)
    skew.ax.axvline(-10, color='darkblue', linestyle='--', alpha=0.7)
    skew.ax.text(-21.5,300, 'Hail Growth\nLayer', rotation=60)
    skew.ax.legend(loc='upper right')
    
    skew.ax.set_title('Hail - Tornado Supercell Average Profile')
    #skew.ax.set_title('Supercell - Non-Supercell Average Profile')
    
    #plt.savefig(r'C:\Users\dylan\Downloads/' + 'SC-NSC-Soundings_Overlay.png',dpi=300)
    #plt.savefig(r'C:\Users\dylan\Downloads/' + 'Hail-Tor-Soundings_Overlay.png',dpi=300)