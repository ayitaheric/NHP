# -*- coding: utf-8 -*-
"""
Created on Fri Jan 19 08:54:50 2024

@author: dylan

PLot centred over AB that gives upper air maps for 250, 500, 700, and 850 hPA, as well as the surface
"""

from netCDF4 import Dataset as NetCDFFile
import metpy.calc as mpcalc
from metpy.units import units
import matplotlib.pyplot as plt
import numpy as np
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from scipy.ndimage import gaussian_filter
from datetime import datetime as dt, timedelta
import os
import pandas as pd

plt.close('all')

#setting up figure
def figure(height):
    fig = plt.figure(figsize = (15,15))
    ax = fig.add_subplot(1, 1, 1, projection=ccrs.Mercator())
    ax.add_feature(cfeature.STATES.with_scale('50m'), linestyle = '-', alpha = 0.3, zorder = 2)
    
    if height == 1000:
        lonW = -122
        lonE = -105
        latS = 60
        latN = 47
    else:
        lonW = -140
        lonE = -95
        latS = 60
        latN = 42
        
    ax.set_extent([lonW, lonE, latS, latN])
    
    return(fig, ax)

def units_func(nc, height, sig_val): 
    lats = nc.variables['latitude'][:] 
    lons = nc.variables['longitude'][:]
    
    level = np.array(nc['level']).tolist().index(height)

    u = units.Quantity(nc['u'][0,level,:,:], "m/s")
    v = units.Quantity(nc['v'][0,level,:,:], "m/s")
    gh = gaussian_filter(nc['z'][0,level,:,:] / 9.80665, sig_val)
    t = gaussian_filter(units.Quantity(nc['t'][0,level,:,:], 'kelvin').to(units('degC')), sig_val)
    w = gaussian_filter(mpcalc.wind_speed(u.to(units.knots),v.to(units.knots)), sig_val)
    
    #creating lats and lons into a gridded format
    lons,lats= np.meshgrid(lons,lats)
    pres = units.Quantity(nc['level'][level], 'hPa')
    q = units.Quantity(nc['q'][0,level,:,:], "kg/kg")
    dew = gaussian_filter(mpcalc.dewpoint_from_specific_humidity(pres,units.Quantity(t, 'degC'), q), sig_val)

    return (lats, lons, level, u, v, w, gh, t, dew)

#instructions for what to do with each upper air map level
def _250mb():
    clevsgh=np.arange(9600,11521,120)
    clevsw=np.arange(20,141, 20)
    cmap = 'BuPu'
    
    return clevsw, clevsgh, cmap

def _300mb():
    clevsgh=np.arange(8280,10440,120)
    clevsw=np.arange(20,141, 20)
    cmap = 'BuPu'
    
    return clevsw, clevsgh, cmap

def _500mb():
    clevsgh=np.arange(5280,6180,60)
    clevsw=np.arange(20,101, 20)
    cmap = 'Blues'
    clevst = np.arange(-45,20, 5)
    
    return clevsw, clevsgh, cmap, clevst

def _700mb():
    clevsgh=np.arange(2100,4020,30)
    clevsw=np.arange(20,60, 10)
    cmap = 'YlOrBr'
    clevsd = np.arange(-4,20, 4)
    clevst = np.arange(-20,32, 4)
    
    return clevsw, clevsgh, cmap, clevsd, clevst

def _850mb():
    clevsgh=np.arange(1050,2010,30)
    clevsw=np.arange(20,60, 10)
    cmap = 'YlOrRd'
    clevsd = np.arange(6,32, 2)
    clevst = np.arange(-4,32, 4)
    
    return clevsw, clevsgh, cmap, clevsd, clevst

def sfc():
    clevsp = np.arange(972,1032,2)
    clevsd = np.arange(-16,32,2)
    clevst = np.arange(-4,40,4)
    clevsm = np.arange(952,1064,4)
    cmap = 'Greens'
    
    return clevsp, clevsd, clevst, clevsm, cmap


filedir = r'C:\Users\dylan\Downloads\NTP Related Documents\2023-2024 Research Projects\NHP Hail Work\ERA5 Data/'
plotsdir = r'C:\Users\dylan\Downloads\NTP Related Documents\2023-2024 Research Projects\NHP Hail Work\Plots/'

directory = os.fsencode(filedir)

df = pd.read_csv(r'C:\Users\dylan\Downloads\NTP Related Documents\2023-2024 Research Projects\NHP Hail Work\Hail_data_MATEUSZ\Chosen Profiles.csv')
df['datetime'] = (pd.to_datetime(df['ERA5_DATE'])).dt.strftime('%Y-%m-%d-%H')

for file in os.listdir(directory):
    
    f = os.fsdecode(file)
    
    date = f[:len(f)-6]
    hour = f[len(f)-5:len(f)-3]

    path = os.path.join(plotsdir, date)
    
    if not os.path.isdir(plotsdir + f'{date}'):
        os.mkdir(path)

    sig_val = 3
    
    #opening NetCDF file and confirming it has all hours (should be length of 1460)
    nc = NetCDFFile(filedir+f)
    
    heights = [250, 500, 700, 850,1000]
    #heights = [1000]
    
    df2 = df[df.datetime == (file[:len(f)-3]).decode("utf-8")]
    df2 = df2.reset_index(drop=True)

    if len(df2) >=1:
        pass

    #iterating through each pressure level for each time stamp
    for height in heights:
        
        lats, lons, level, u, v, wind, gh, temp, dew = units_func(nc, height, sig_val)
        
        if height == 250:
            clevsw, clevsgh, cmap = _250mb()
        elif height == 300:
            clevsw, clevsgh, cmap,  clevst = _300mb()
        elif height == 500:
            clevsw, clevsgh, cmap,  clevst = _500mb()
        elif height == 700:
            clevsw, clevsgh, cmap, clevsd, clevst = _700mb()
        elif height == 850:
            clevsw, clevsgh, cmap, clevsd, clevst = _850mb()
        elif height == 1000:
            clevsp, clevsd, clevst, clevsm, cmap = sfc()
        
        #plotting SLP contour lines and the low pressure
        if height != 1000:
            fig, ax = figure(height)
            
            c = ax.contourf(lons,lats,wind,clevsw, zorder=0, cmap=cmap, transform=ccrs.PlateCarree())
            cb = plt.colorbar(c, ax=ax, orientation='horizontal', pad=0.01)
            cb.set_label('Wind Speed (knots)', fontsize=15, labelpad=10)
            
            c = ax.contour(lons,lats,gh,clevsgh, zorder=3, colors = 'black', transform=ccrs.PlateCarree())
            plt.clabel(c, fmt = '%1.0f', inline_spacing=14)
            time = dt(1900,1,1) + timedelta(hours = int(nc["time"][0]))
            plt.title(f'{time.strftime("%Y-%m-%d %H UTC")} - {nc["level"][level]} hPa', fontsize = 15)
        
        if height in (700,850):
            c = ax.contour(lons,lats,temp,clevst, zorder=3, colors = 'red', transform=ccrs.PlateCarree())
            plt.clabel(c, fmt = '%1.0f', inline_spacing=14)
            if height in (700, 850):
                c = ax.contour(lons,lats,dew,clevsd, zorder=3, colors = 'green', transform=ccrs.PlateCarree())
                plt.clabel(c, fmt = '%1.0f', inline_spacing=14)
                
        if height == 1000:
            fig, ax = figure(height)
            
            u_nc = NetCDFFile(r'C:\Users\dylan\Downloads\NTP Related Documents\2023-2024 Research Projects\NHP Hail Work\ERA5 Complete Data/' + f'ERA5-{date}-{hour}-u-sfc.nc')
            v_nc = NetCDFFile(r'C:\Users\dylan\Downloads\NTP Related Documents\2023-2024 Research Projects\NHP Hail Work\ERA5 Complete Data/' + f'ERA5-{date}-{hour}-v-sfc.nc')
            d_nc = NetCDFFile(r'C:\Users\dylan\Downloads\NTP Related Documents\2023-2024 Research Projects\NHP Hail Work\ERA5 Complete Data/' + f'ERA5-{date}-{hour}-dew-sfc.nc')
            m_nc = NetCDFFile(r'C:\Users\dylan\Downloads\NTP Related Documents\2023-2024 Research Projects\NHP Hail Work\ERA5 Complete Data/' + f'ERA5-{date}-{hour}-mslp.nc')

            lats = u_nc['latitude'][:]
            lons = u_nc['longitude'][:]
            lons,lats= np.meshgrid(lons,lats)
            
            u = units.Quantity(u_nc['u10'][0,:,:], 'm/s').to(units['km/h'])
            v = units.Quantity(v_nc['v10'][0,:,:], 'm/s').to(units['km/h'])
            d = gaussian_filter(units.Quantity(d_nc['d2m'][0,:,:], 'kelvin').to(units('degC')), sig_val)
            m = gaussian_filter(units.Quantity(m_nc['msl'][0,:,:], 'Pa').to(units['hPa']), sig_val)
            #wind = mpcalc.wind_speed(u.to(units.knots),v.to(units.knots))
            
            c = ax.contour(lons,lats,d,clevsd, zorder=0, cmap=cmap, transform=ccrs.PlateCarree())
            plt.clabel(c, fmt = '%1.0f', inline_spacing=14)
            
            c = ax.contour(lons,lats,m,clevsm, zorder=0, colors='black', transform=ccrs.PlateCarree())
            plt.clabel(c, fmt = '%1.0f', inline_spacing=14)
            
            ax.barbs(lons[::7,::7], lats[::7,::7], u.magnitude[::7,::7], v.magnitude[::7,::7], transform=ccrs.PlateCarree())
            
        ax.plot(df2.Longitude[0], df2.Latitude[0], '*', markersize=10, c='black', transform=ccrs.PlateCarree(), zorder=2)
        
        plt.title(f'{date} {hour} UTC - {height} hPa Map')
        #plt.savefig(plotsdir + f'{date}/' + f'{height}mb_UA_Map-{hour}Z.png')
