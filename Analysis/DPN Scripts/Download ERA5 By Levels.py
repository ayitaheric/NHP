# -*- coding: utf-8 -*-
"""
Created on Sun Jul 16 14:31:32 2023

@author: wangm/dylan

This script follows the Copernicus Climate Data Store instruction on how to download the ERA5 data.
First, need to register, get an API access key, and install the CDS API client.

"""

import cdsapi
import pandas as pd
import numpy as np
import sys, os

def Retrieve_Data(datetime, project_dir, folder_dir):
    year = datetime.strftime("%Y")
    month = datetime.strftime("%m")
    day = datetime.strftime("%d")
    hour = datetime.strftime("%H")
    
    c.retrieve(
        'reanalysis-era5-pressure-levels',       # this is the name of the data on the CDS store
        {
            'product_type': 'reanalysis',
            'variable': [
                'divergence', 'fraction_of_cloud_cover', 'geopotential',
                'potential_vorticity', 'relative_humidity', 
                'specific_humidity', 'temperature','u_component_of_wind', 
                'v_component_of_wind','vertical_velocity','vorticity',
                        ],
        'pressure_level': ['100','125',
                           '150', '175', '200',
                           '225', '250', '300',
                           '350', '400', '450',
                           '500', '550', '600',
                           '650', '700', '750',
                           '775', '800', '825',
                           '850', '875', '900',
                           '925', '950', '975',
                           '1000',],
        'year': f'{year}',
        'month': f'{month}',
        'day': f'{day}',
        'time': f'{hour}' + ':00',
        'area': [65,-145, 40,-90,],
        'format': 'netcdf',
       },
        project_dir + folder_dir + f'{year}-{month}-{day}-{hour}'+'.nc')


# Read in individual files that contain the date and time of each event wanted.
# Here we just want to extract the date and time of the event and feed that into the API
# request.

project_dir = r'C:\Users\dylan\Downloads\NHP Hail Work/'
folder_dir = 'ERA5 Data/'
        
data_event = pd.read_csv(project_dir + 'Hail_data_MATEUSZ\Chosen Profiles.csv')

data_event['datetime'] = pd.to_datetime(data_event['ERA5_DATE'])

event_lat = np.array(data_event.Latitude,dtype='float')
event_lon = np.array(data_event.Longitude,dtype='float')

# initialize the API client

c = cdsapi.Client()

# specify a directory to save the data to

# loop through each event listed in the CSV file. The CDS API will retrieve and save the ERA5 reanalysis
# data for each event.

for f in range(0,len(event_lat)):
    file_name = data_event.datetime[f].strftime('%Y-%m-%d-%H')

    if os.path.isfile(project_dir + folder_dir + f'{file_name}.nc'):
        pass
    else:
        Retrieve_Data(data_event.datetime[f], project_dir, folder_dir)
