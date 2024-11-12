# -*- coding: utf-8 -*-
"""
Created on Tue Jan 16 12:31:43 2024

@author: dylan

Map showing where supercell and non-supercell hail events were recorded. This is overlayed over
Alberta's elevation in m (ASL).
"""
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import matplotlib.ticker as mticker
from datetime import datetime as dt
import matplotlib.colors as colors
from netCDF4 import Dataset as NetCDFFile
from cartopy.mpl.ticker import (LongitudeFormatter, LatitudeFormatter,
                                LatitudeLocator, LongitudeLocator)

def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=100):
    new_cmap = colors.LinearSegmentedColormap.from_list(
        'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval),
        cmap(np.linspace(minval, maxval, n)))
    return new_cmap

df = pd.read_csv(r'C:\Users\dylan\Downloads\NTP Related Documents\2023-2024 Research Projects\NHP Hail Work\Hail_data_MATEUSZ\Chosen Profiles.csv')

plt.close('all')

fig = plt.figure(figsize = (10,8))
ax = fig.add_subplot(1, 1, 1, projection=ccrs.Mercator())

#adding background cartopy features
ax.add_feature(cfeature.STATES.with_scale('50m'), linestyle = '-', alpha = 0.5, zorder = 2)

ax.set_xticks([-120,-115,-110], crs=ccrs.PlateCarree())
ax.set_yticks([50,52,54], crs=ccrs.PlateCarree())
lon_formatter = LongitudeFormatter(zero_direction_label=True)
lat_formatter = LatitudeFormatter()
ax.xaxis.set_major_formatter(lon_formatter)
ax.yaxis.set_major_formatter(lat_formatter)

#map boundaries, which covers all of Canada
longW = -120.5
longE = -109.5
latS = 48.9
latN = 55

ax.set_extent([longW, longE, latS, latN])

#finding x and y values of supercell and non-supercell events
xsc = df.loc[df['Storm Type'] == 'Y', 'Longitude']
ysc = df.loc[df['Storm Type'] == 'Y', 'Latitude']

xnsc = df.loc[df['Storm Type'] == 'N', 'Longitude']
ynsc = df.loc[df['Storm Type'] == 'N', 'Latitude']

#plotting events
ax.scatter(xsc, ysc, marker = 'o', color = 'red', zorder = 1, s = 30, transform = ccrs.PlateCarree(), edgecolor='black', label = 'Supercell')
ax.scatter(xnsc, ynsc, marker = 'o', color = 'blue', zorder = 1, s = 30, transform = ccrs.PlateCarree(), edgecolor='black', label = 'Non-Supercell')

#plotting of cities and Canada's record hail stone
ax.scatter(-114.2874, 52.1809, marker = '*', color = 'red', zorder = 1, s = 105, transform = ccrs.PlateCarree(), edgecolor='black', label = 'Record Hailstone')
ax.scatter(-113.4937, 53.5461, marker = 's', color = 'black', zorder = 1, s = 50, transform = ccrs.PlateCarree(), edgecolor='black', label = 'Edmonton')
ax.scatter(-114.1199, 51.7951, marker = 's', color = 'cyan', zorder = 1, s = 50, transform = ccrs.PlateCarree(), edgecolor='black', label = 'Olds')
ax.scatter(-114.0719, 51.0447, marker = 's', color = 'grey', zorder = 1, s = 50, transform = ccrs.PlateCarree(), edgecolor='black', label = 'Calgary')

plt.legend(loc="upper right", framealpha = 1, prop={'size': 12})
plt.title('NHP Hail Sample Locations', size = 15, weight = 'bold')

#file directory and specific file that will be examined
filedir = r'C:\Users\Dylan\Documents\Research Projects\Historical Winter Storms\Various files/'
file = 'Terrain_Height.nc'

#opening NetCDF file and confirming it has all hours (should be length of 1460)
nc = NetCDFFile(filedir+file)

lat = (nc.variables['XLAT'][:][:].squeeze())[700:1500,200:600]
lon = (nc.variables['XLONG'][:][:].squeeze())[700:1500,200:600]
hgt = (nc.variables['HGT'][:][:].squeeze())[700:1500,200:600]

#adding background cartopy features
ax.add_feature(cfeature.STATES.with_scale('50m'), linestyle = '-', alpha = 0.5, zorder = 2)

clevs=np.arange(0,3001, 100)

#colour map for terrain
cmap = plt.cm.terrain

new_cmap = truncate_colormap(cmap, 0.25, 1)

c = ax.contourf(lon,lat,hgt,clevs, zorder=0, cmap = new_cmap, transform=ccrs.PlateCarree())

#making a colourbar label on the side of the figure for dBZ
cbar_ax = fig.add_axes([0.90, 0.11, 0.02, 0.771])
cbar = plt.colorbar(c, cax = cbar_ax)
cbar.ax.set_ylabel('Terrain Height (m)', weight = 'bold', size = 14)
cbar.set_ticks(np.arange(0, 3001, 500))
cbar.ax.set_yticklabels(np.arange(0, 3001, 500), weight = 'bold', size=14)

#plt.savefig(r'C:\Users\dylan\Downloads\NHP Hail Work\Plots\NHP Hail Report Locations.png', bbox_inches='tight', dpi=300)