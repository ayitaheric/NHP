# -*- coding: utf-8 -*-
"""
Created on Tue Feb  6 12:21:41 2024

@author: dylan

Script plotting two skew-ts side by side. Features wind barbs as well as the ML LCL, LFC, EL, and CAPE
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

#finding pressure intervals
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

warnings.filterwarnings("ignore", category=DeprecationWarning) 
warnings.filterwarnings("ignore", category=FutureWarning)

plt.close('all')

folder_dir = r'C:\Users\dylan\Downloads\NTP Related Documents\2023-2024 Research Projects\NHP Hail Work\SC_unrotated/'
#arr = np.empty((5,), dtype=np.ndarray)

fig = plt.figure(figsize=(12,6))
skew = [np.nan]*2
skew[0] = SkewT(fig, subplot = (1,2,1))
skew[1] = SkewT(fig, subplot = (1,2,2))

#going over both datasets
for h in range(2):
    if h ==0:
        df = pd.read_csv(folder_dir + 'SC_RAW_data.csv')
        title = 'Supercell'
        #title = 'Hail'
    elif h ==1:
        df = pd.read_csv(folder_dir + 'NSC_RAW_data.csv')
        title = 'Non-Supercell'
        #df = pd.read_csv(folder_dir + 'Region1_AB_RAW_data.csv')
        #title = 'Tornado'
        
    a = df.altitude
    p = df.pressure
    t = df.temp2
    d = df.dpt2
    
    #finding u and v components for calculations
    u = []
    v = []
    data_tv = []
    for i in range(len(df.ws)):
        u.append((mpcalc.wind_components(df.ws[i] * units('knots').to(units('m/s')),df.wd[i]*units.deg)[0]).magnitude)
        v.append((mpcalc.wind_components(df.ws[i] * units('knots').to(units('m/s')),df.wd[i]*units.deg)[1]).magnitude)
        data_tv.append(thermo.virtemp(p[i], t[i], d[i]))
        
    data_tv = np.array(data_tv)
    a = [x - a[0] for x in a]

    wdir = mpcalc.wind_direction(u*units('m/s'),v*units('m/s'))
    wspd = mpcalc.wind_speed(u*units('m/s'),v*units('m/s'))
    
    data_profset = profile.create_profile(profile='default', pres=p, hght=a,
                                          tmpc=t, dwpc=d, wdir=wdir, wspd=wspd)

    sfcpres = data_profset.pres[data_profset.sfc]
    # sfcpres = profset.pres[profset.sfc]

    # p3km = interp.pres(prof, interp.to_msl(prof, 3000.))
    
    #parcel calculations
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
    
    #plotting data on skew-t
    skew[h].plot(p,t,c=[1,0,0],linewidth=1,ls='-')
    skew[h].plot(p,data_tv, c=[1,0,0], linewidth=1,ls=':')
    skew[h].plot(p,d,c=[0,1,0],linewidth=1,ls='-')
    # skews.plot(data_sfcpcl.ptrace,data_sfcpcl.ttrace,color=[0.5,0,0.5],linewidth=1, linestyle='--')
    skew[h].plot(data_mlpcl.ptrace,data_mlpcl.ttrace,color=[1,0,1],linewidth=1, linestyle='-')
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
    
    #plotting ML LCL, LFC, EL
    skew[h].plot([lcl_pres, lcl_pres], [lcl_temp-1, lcl_temp+1], 'k-', lw=1)
    skew[h].ax.text(lcl_temp.data+2, lcl_pres, 'ML LCL',fontsize=8)
    skew[h].plot([lfc_pres, lfc_pres], [lfc_temp-1, lfc_temp+1], 'k-', lw=1)
    skew[h].ax.text(lfc_temp.data+2, lfc_pres, 'ML LFC',fontsize=8)
    skew[h].plot([el_pres, el_pres], [el_temp-1, el_temp+1], 'k-', lw=1)
    skew[h].ax.text(el_temp.data+2, el_pres, 'ML EL',fontsize=8)
    
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
    skew[h].shade_cin(data_pres_cinshade, data_temp_cinshade, parc_tv_cinshade)
    skew[h].shade_cape(p, data_tv, parc_tvshade)
    
    #p_,u_,v_ = pressure_interval(p.magnitude,u,v)
    p_,u_,v_ = pressure_interval(p,u,v)
    skew[h].plot_barbs(p_,u_,v_)

    # skews[fin].shade_cin(env_pshade.data, env_tvshade.data, data_mlpcl.ttrace.data)
    # skews[fin].shade_cape(env_pshade, env_tvshade, data_mlpcl.ttrace.data)
    
    # sys.exit()

    skew[h].ax.set_xlim(-40,50)
    skew[h].ax.set_ylim(1000,100)
    skew[h].ax.set_xticks(np.arange(-40,50.01,10))
    skew[h].ax.set_yticks(np.arange(1000,99.99,-100))
    skew[h].ax.set_xticklabels(['-40','-30','-20','-10','0','10','20','30','40','50'],fontsize=11)
    if h == 0:
        skew[h].ax.set_yticklabels(['1000','900','800','700','600','500','400','300','200','100'],fontsize=11)
    else:
        skew[h].ax.set_yticklabels([],fontsize=11)
    skew[h].ax.set_xlabel('Temperature ($\mathregular{^o}$C)',fontsize=12)
    if h == 0:
        skew[h].ax.set_ylabel('Pressure (hPa)', fontsize=12)
    #skew[h].ax.axvline(0,color='c',linestyle='--',linewidth=1.5)
    skew[h].plot_dry_adiabats(t0=(np.arange(-50,260,10))*units('degC'),color='gray',linewidth=0.5,linestyle='--', alpha = 0.3)
    skew[h].plot_moist_adiabats(color='gray',linewidth=0.5,linestyle='--', alpha = 0.3)
    skew[h].plot_mixing_lines(color='gray',linewidth=0.5,linestyle=':')
    skew[h].ax.axvline(-30, color='darkblue', linestyle='--', alpha=0.7)
    skew[h].ax.axvline(-10, color='darkblue', linestyle='--', alpha=0.7)
    skew[h].ax.text(-21.5,300, 'Hail Growth\nLayer', rotation=60)
    skew[h].ax.set_title(f'{title}')
    
    #plt.savefig(r'C:\Users\dylan\Downloads\NHP Hail Work/' + 'SC-NSC Separate Skew-ts.png', dpi=300)