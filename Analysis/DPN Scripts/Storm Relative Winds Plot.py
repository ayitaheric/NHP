# -*- coding: utf-8 -*-
"""
Created on Wed Mar  6 22:01:34 2024

@author: dylan

Plot showing BRM for two different datasets (in this case supercell and non-supercell) of the 0-6km 
mean wind storm relative winds. There are 9 different plots, with each row consisting of u and v 
components, and thethird of the U component. The coloumns are split into SC, NSC, and both.
"""
import metpy.calc as mpcalc
from metpy.units import units, pandas_dataframe_to_unit_arrays
import pandas as pd
import matplotlib.pyplot as plt
from statistics import mean, quantiles
import numpy as np
import os

plt.close('all')

plot_dir = r'C:\Users\dylan\Downloads\NTP Related Documents\2023-2024 Research Projects\NHP Hail Work\Plots/'
folder_dir = r'C:\Users\dylan\Downloads\NTP Related Documents\2023-2024 Research Projects\NHP Hail Work\Hail_data_MATEUSZ\Raw_profiles/'

df_og = pd.read_csv(r'C:\Users\dylan\Downloads\NTP Related Documents\2023-2024 Research Projects\NHP Hail Work\Hail_data_MATEUSZ\Chosen Profiles.csv')
df_og['datetime'] = pd.to_datetime(df_og['ERA5_DATE'])

fig, ax = plt.subplots(3, 3, figsize=(18,12))
plt.setp(ax, ylim=(0,12000))

#iterating over both storm types
for k in range(2):
    if k == 0:
        df = df_og[df_og['Storm Type'] == 'Y']
        c = 'red'
    elif k == 1:
        df = df_og[df_og['Storm Type'] == 'N']
        c = 'blue'
        
    df = df.reset_index()
    
    sr_us = []
    sr_vs = []
    wind_speeds = []

    for i in range(len(df['Event Number'])):      
        time = df.datetime[i].strftime('_%Y-%m-%d_')
        
        if int(df['Event Number'][i]) <10:
            num = f'_000{df["Event Number"][i]}'
        else:
            num = f'_00{df["Event Number"][i]}'
            
        hour = f'{df.datetime[i].strftime("%H")}_'
        
        df_prof = pd.read_csv(folder_dir + f'HAIL{num}{time}{hour}UTC.csv')
    	#pressure	altitude	temp	dpt	wd	ws	u	v	Lat	Lon	Date

        unit_types = {'pressure':'hPa', 'temp':'degC', 'dpt':'degC', 'wd':'deg', 'ws':'m/s',\
                      'u':'m/s', 'v':'m/s'}
        united_array=pandas_dataframe_to_unit_arrays(df_prof,column_units=unit_types)

        a = united_array['altitude']
        p = united_array['pressure']
        u = united_array['u']
        v = united_array['v']
        a = [x - a[0] for x in a]
        
        #bunkers right mover
        if k == 0:
            brm = (mpcalc.bunkers_storm_motion(p, u, v, height = a*units.m)[0]).magnitude
        elif k == 1:
            #brm = (mpcalc.bunkers_storm_motion(p, u, v, height = a*units.m)[0]).magnitude
            
            brm1 = (mpcalc.bunkers_storm_motion(p, u, v, height = a*units.m)[0]).magnitude
            brm2 = (mpcalc.bunkers_storm_motion(p, u, v, height = a*units.m)[2]).magnitude

            brm[0] = (brm1[0]+brm2[0])/2
            brm[1] = (brm1[1]+brm2[1])/2
            

        #BRM wind direction
        wd = mpcalc.wind_direction(brm[0]*units('m/s'),brm[1]*units('m/s'))
        #print(wd.magnitude)
        #wind difference between BRM and 270 deg (x-axis)
        wd_diff = 270 - wd.magnitude
        #wind dir and speed of profile
        wd_prof = mpcalc.wind_direction(u,v)
        ws_prof = mpcalc.wind_speed(u,v)
        #new wind dir when rotating hodograph
        wd_new = wd_prof.magnitude + wd_diff
        #print(wd_new)
        #new BRM values
        brm_ws = mpcalc.wind_speed(brm[0]*units('m/s'),brm[1]*units('m/s'))
        #print(brm_ws.magnitude)
        brm_u,brm_v = mpcalc.wind_components(brm_ws,270*units.deg)

        #new u and v components
        u_new,v_new = mpcalc.wind_components(ws_prof,wd_new*units.deg)

        sr_us.append(u_new.magnitude - brm_u.magnitude)
        sr_vs.append(v_new.magnitude - brm_v.magnitude)
        
        wind_speeds.append((mpcalc.wind_speed(((df_prof.u - brm[0]).to_list())*units('m/s'),((df_prof.v - brm[1]).to_list())*units('m/s')))/units('m/s'))
        
    sr_u = (list(map(mean, zip(*sr_us))))
    sr_v = (list(map(mean, zip(*sr_vs))))
    wind_speed = (list(map(mean, zip(*wind_speeds))))
    sr_u_min = np.array(list(map(quantiles, zip(*sr_us))))
    sr_u_max = np.array(list(map(quantiles, zip(*sr_us))))
    sr_v_min = np.array(list(map(quantiles, zip(*sr_vs))))
    sr_v_max = np.array(list(map(quantiles, zip(*sr_vs))))
    wind_speed_min = np.array(list(map(quantiles, zip(*wind_speeds))))
    wind_speed_max = np.array(list(map(quantiles, zip(*wind_speeds))))
    
    ax[0,k].plot(sr_u,a,color=c)
    ax[0,2].plot(sr_u,a,color=c)
    ax[1,k].plot(sr_v,a,color=c)
    ax[1,2].plot(sr_v,a,color=c)
    ax[0,k].fill_betweenx(a, sr_u_min[:,0], sr_u_max[:,2], facecolor=c, alpha=0.2)
    ax[1,k].fill_betweenx(a, sr_v_min[:,0], sr_v_max[:,2], facecolor=c, alpha=0.2)
    ax[0,2].fill_betweenx(a, sr_u_min[:,0], sr_u_max[:,2], facecolor=c, alpha=0.2)
    ax[1,2].fill_betweenx(a, sr_v_min[:,0], sr_v_max[:,2], facecolor=c, alpha=0.2)

    ax[2,k].plot(wind_speed,a,color=c)
    ax[2,k].fill_betweenx(a, wind_speed_min[:,0], wind_speed_max[:,2], facecolor=c, alpha=0.2)
    ax[2,2].plot(wind_speed,a,color=c)
    ax[2,2].fill_betweenx(a, wind_speed_min[:,0], wind_speed_max[:,2], facecolor=c, alpha=0.2)
    
    ax[k,0].grid()
    ax[k,1].grid()
    ax[k,2].grid()
    ax[k,0].set_ylabel('Height (m AGL)')
    ax[k,1].set_ylabel('Height (m AGL)')
    ax[k,2].set_ylabel('Height (m AGL)')
    ax[k,0].set_xlabel('Speed (m s\N{SUPERSCRIPT MINUS}\N{SUPERSCRIPT ONE})')
    ax[k,1].set_xlabel('Speed (m s\N{SUPERSCRIPT MINUS}\N{SUPERSCRIPT ONE})')
    ax[k,2].set_xlabel('Speed (m s\N{SUPERSCRIPT MINUS}\N{SUPERSCRIPT ONE})')
    #ax[k,0].set_xlim(-25,45)
    #ax[k,1].set_xlim(-25,45)
    #ax[k,2].set_xlim(-25,45)
    if k == 0:
        ax[k,0].set_xlim(-20,25)
        ax[k,1].set_xlim(-20,25)
        ax[k,2].set_xlim(-20,25)
    else:
        ax[k,0].set_xlim(-10,20)
        ax[k,1].set_xlim(-10,20)
        ax[k,2].set_xlim(-10,20)
    
"""
for k in range(1):
    k=1
    if k == 1:
        c = 'blue'
    
    sr_us = []
    sr_vs = []
    wind_speeds = []
"""   
    #folder_dir = r'C:\Users\dylan\Downloads\NHP Hail Work\TORNADO_raw_profiles\AB/'
"""  
    directory = os.fsencode(folder_dir)
    
    for file in os.listdir(directory):
        f = os.fsdecode(file)
        
        df_prof = pd.read_csv(folder_dir + f)       
    	#pressure	altitude	temp	dpt	wd	ws	u	v	Lat	Lon	Date

        unit_types = {'temp':'degC','dpt':'degC','pres':'hPa', 'wd':'deg', 'ws':'m/s'}
        united_array=pandas_dataframe_to_unit_arrays(df_prof,column_units=unit_types)

        a = united_array['alt']
        p = united_array['pres']
        ws = united_array['ws']
        wd = united_array['wd']
        a = [x - a[0] for x in a]
        
        u,v = mpcalc.wind_components(ws,wd)
        
        #bunkers right mover
        brm = (mpcalc.bunkers_storm_motion(p, u, v, height = a*units.m)[0]).magnitude
        
        #BRM wind direction
        wd = mpcalc.wind_direction(brm[0]*units('m/s'),brm[1]*units('m/s'))
        #wind difference between BRM and 270 deg (x-axis)
        wd_diff = 270 - wd.magnitude
        #wind dir and speed of profile
        wd_prof = mpcalc.wind_direction(u,v)
        ws_prof = mpcalc.wind_speed(u,v)
        #new wind dir when rotating hodograph
        wd_new = wd_prof.magnitude + wd_diff
        #print(wd_new)
        #new BRM values
        brm_ws = mpcalc.wind_speed(brm[0]*units('m/s'),brm[1]*units('m/s'))
        brm_u,brm_v = mpcalc.wind_components(brm_ws,270*units.deg)

        #new u and v components
        u_new,v_new = mpcalc.wind_components(ws_prof,wd_new*units.deg)

        sr_us.append((u_new.magnitude - brm_u.magnitude))
        sr_vs.append((v_new.magnitude - brm_v.magnitude))
        
        wind_speeds.append((mpcalc.wind_speed(((u.magnitude - brm[0]))*units('m/s'),((v.magnitude - brm[1]))*units('m/s')))/units('m/s'))
        
    sr_u = (list(map(mean, zip(*sr_us))))
    sr_v = (list(map(mean, zip(*sr_vs))))
    wind_speed = (list(map(mean, zip(*wind_speeds))))
    sr_u_min = np.array(list(map(quantiles, zip(*sr_us))))
    sr_u_max = np.array(list(map(quantiles, zip(*sr_us))))
    sr_v_min = np.array(list(map(quantiles, zip(*sr_vs))))
    sr_v_max = np.array(list(map(quantiles, zip(*sr_vs))))
    wind_speed_min = np.array(list(map(quantiles, zip(*wind_speeds))))
    wind_speed_max = np.array(list(map(quantiles, zip(*wind_speeds))))
    
    ax[0,k].plot(sr_u,a,color=c)
    ax[0,2].plot(sr_u,a,color=c)
    ax[1,k].plot(sr_v,a,color=c)
    ax[1,2].plot(sr_v,a,color=c)
    ax[0,k].fill_betweenx(a, sr_u_min[:,0], sr_u_max[:,2], facecolor=c, alpha=0.2)
    ax[1,k].fill_betweenx(a, sr_v_min[:,0], sr_v_max[:,2], facecolor=c, alpha=0.2)
    ax[0,2].fill_betweenx(a, sr_u_min[:,0], sr_u_max[:,2], facecolor=c, alpha=0.2)
    ax[1,2].fill_betweenx(a, sr_v_min[:,0], sr_v_max[:,2], facecolor=c, alpha=0.2)

    ax[2,k].plot(wind_speed,a,color=c)
    ax[2,k].fill_betweenx(a, wind_speed_min[:,0], wind_speed_max[:,2], facecolor=c, alpha=0.2)
    ax[2,2].plot(wind_speed,a,color=c)
    ax[2,2].fill_betweenx(a, wind_speed_min[:,0], wind_speed_max[:,2], facecolor=c, alpha=0.2)
    
    ax[k,0].grid()
    ax[k,1].grid()
    ax[k,2].grid()
    ax[k,0].set_ylabel('Height (m AGL)')
    ax[k,1].set_ylabel('Height (m AGL)')
    ax[k,2].set_ylabel('Height (m AGL)')
    ax[k,0].set_xlabel('Speed (m s\N{SUPERSCRIPT MINUS}\N{SUPERSCRIPT ONE})')
    ax[k,1].set_xlabel('Speed (m s\N{SUPERSCRIPT MINUS}\N{SUPERSCRIPT ONE})')
    ax[k,2].set_xlabel('Speed (m s\N{SUPERSCRIPT MINUS}\N{SUPERSCRIPT ONE})')
    ax[k,0].set_xlim(-20,30)
    ax[k,1].set_xlim(-20,30)
    ax[k,2].set_xlim(-20,30)
"""  
ax[2,0].grid()
ax[2,1].grid()
ax[2,2].grid()
ax[2,0].set_ylabel('Height (m AGL)')
ax[2,1].set_ylabel('Height (m AGL)')
ax[2,2].set_ylabel('Height (m AGL)')
ax[2,0].set_xlabel('Speed (m s\N{SUPERSCRIPT MINUS}\N{SUPERSCRIPT ONE})')
ax[2,1].set_xlabel('Speed (m s\N{SUPERSCRIPT MINUS}\N{SUPERSCRIPT ONE})')
ax[2,2].set_xlabel('Speed (m s\N{SUPERSCRIPT MINUS}\N{SUPERSCRIPT ONE})')
#ax[2,0].set_xlim(0,50)
#ax[2,1].set_xlim(0,50)
#ax[2,2].set_xlim(0,50)
ax[2,0].set_xlim(0,30)
ax[2,1].set_xlim(0,30)
ax[2,2].set_xlim(0,30)
    
ax[0,0].set_title('Supercell',size=12,weight='bold')
ax[0,1].set_title('Non-Supercell',size=12,weight='bold')
ax[0,2].set_title('Both Storm Types',size=12,weight='bold')

#ax[0,0].set_title('Hail',size=12,weight='bold')
#ax[0,1].set_title('Tornado',size=12,weight='bold')
#ax[0,2].set_title('Both Storm Types',size=12,weight='bold')

titles = ['u', 'v', 'U']

num = 0
for a in ax:
    for b in a:
        b.annotate(f'{titles[num]}'+r'$_{_{SR}}$', xy=(0, 1), xycoords='axes fraction', fontsize=16, xytext=(5, -5), textcoords='offset points',
                horizontalalignment='left', verticalalignment='top')
    num+=1
plt.suptitle('SC BRM - NSC BRM/0-6km Mean Wind Storm Relative Winds',y=0.93,fontsize=16,weight='bold') 
#plt.savefig(r'C:\Users\dylan\Downloads\NHP Hail Work\Plots\SC-NSC-SR_Winds_Plot-BRM0-6kmMeanWind.png',dpi=300, bbox_inches='tight')
#plt.suptitle('Hail - Tornado Supercell Storm Relative Winds',y=0.93,fontsize=16,weight='bold') 
#plt.savefig(r'C:\Users\dylan\Downloads\NHP Hail Work\Plots\Hail-Tornado-SR_Winds_Plot.png',dpi=300, bbox_inches='tight')