# -*- coding: utf-8 -*-
"""
Created on Thu Jan 18 08:51:55 2024

@author: dylan
"""
import cdsapi
import cdsapi
import pandas as pd
import numpy as np
import sys, os
from datetime import datetime as dt, timedelta
import datetime

#definition to retrieve data
def Retrieve_Data(rounded, project_dir, folder_dir, variable):
    date = rounded.strftime('%Y-%m-%d')
    hour = int(rounded.strftime('%H'))

    c = cdsapi.Client()
    c.retrieve('reanalysis-era5-complete', { # Requests follow MARS syntax
                                             # Keywords 'expver' and 'class' can be dropped. They are obsolete
                                             # since their values are imposed by 'reanalysis-era5-complete'
        'date'    : f'{date}',            # The hyphens can be omitted
        'levelist': '1/2/3/4/5/6/7/8/9/10/11/12/13/14/15/16/17/18/19/20/21/22/23/24/25/26/27/28/29/30/31/32/33/34/35/36/37/38/39/40/41/42/43/44/45/46/47/48/49/50/51/52/53/54/55/56/57/58/59/60/61/62/63/64/65/66/67/68/69/70/71/72/73/74/75/76/77/78/79/80/81/82/83/84/85/86/87/88/89/90/91/92/93/94/95/96/97/98/99/100/101/102/103/104/105/106/107/108/109/110/111/112/113/114/115/116/117/118/119/120/121/122/123/124/125/126/127/128/129/130/131/132/133/134/135/136/137',          # 1 is top level, 137 the lowest model level in ERA5. Use '/' to separate values.
        #"expver": "1", #not usually here for non surface
        'levtype' : 'ml',#ml #sfc
        'param'   : '54', # Full information at https://apps.ecmwf.int/codes/grib/param-db/
                                             # The native representation for temperature is spherical harmonics
        'stream'  : 'oper',                  # Denotes ERA5. Ensemble members are selected by 'enda'
        #'00/to/23/by/6', 
        'time'    : f'{hour}',         # You can drop :00:00 and use MARS short-hand notation, instead of '00/06/12/18'
        'type'    : 'an',
        'area'    : '60/-122/47/-105',          # North, West, South, East. Default: global
        'grid'    : '0.25/0.25',               # Latitude/longitude. Default: spherical harmonics or reduced Gaussian grid
        'format'  : 'netcdf',                # Output needs to be regular lat-lon, so only works in combination with 'grid'!
        }, project_dir + folder_dir + f'ERA5-{date}-{hour}-{variable}.nc')
        #}, r'C:\Users\dylan\Downloads\NHP Hail Work\ERA5 Complete Data/'+f'ERA5-{date}-{hour}-{variable}.nc') 


project_dir = r'C:\Users\dylan\Downloads\NHP Hail Work/'
folder_dir = 'ERA5 Complete Data/'

data_event = pd.read_csv(project_dir + 'Hail_data_MATEUSZ\Chosen Profiles.csv')

#variable you are looking to retrieve from data (can have multiple in one netcdf file if needed)
variable = 'p'

#date from csv file iterating through
data_event['datetime'] = pd.to_datetime(data_event['ERA5_DATE'])

#lat and lon of event
event_lat = np.array(data_event.Latitude,dtype='float')
event_lon = np.array(data_event.Longitude,dtype='float')

# initialize the API client

c = cdsapi.Client()

# specify a directory to save the data to

# loop through each event listed in the CSV file. The CDS API will retrieve and save the ERA5 reanalysis
# data for each event.

for f in range(0,len(event_lat)):
    file_name = data_event.datetime[f].strftime('%Y-%m-%d-%H')

    if os.path.isfile(project_dir + folder_dir + f'ERA5-{file_name}-{variable}.grib'):
        pass
    else:
        Retrieve_Data(data_event.datetime[f], project_dir, folder_dir, variable)