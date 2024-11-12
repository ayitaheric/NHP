# Copyright 2016 Cambridge Environmental Research Consultants Ltd.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0
#
# **************************************************************************
# Function      : compute_geopotential_on_ml_netcdf
#
# Author (date) : Mark Jackson (1/12/2016)
#
# Category      : COMPUTATION
#
# OneLineDesc   : Computes geopotential and height on model levels using netCDF files
#
# Description   : Computes geopotential on model levels using netCDF files.
#                 Based on the Python script by Cristian Simarro which uses GRIB files:
#                   https://confluence.ecmwf.int/display/GRIB/Compute+geopotential+on+model+levels
#                 Which was based on the code of the Metview function mvl_geopotential_on_ml:
#                   https://confluence.ecmwf.int/metview/mvl_geopotential_on_ml  
#                 This in turn was based on code from Nils Wedi, the IFS documentation:
#                   https://confluence.ecmwf.int/display/IFS/CY41R1+Official+IFS+Documentation
#                 part III. Dynamics and numerical procedures
#                 and an optimised implementation by Dominique Lucas.
#                 Ported to Python by Cristian Simarro
#
# Parameters    : FileA.nc  - netCDF file with the levelist of t and q. Does not require all levels 
#                             but does require a contiguous set of levels all the way to the bottom.
#                 FileB.nc  - netCDF file with levelist 1 for params z and lnsp
#
# Return Value  : outputs CSV files
#                 with geopotential and height (relative to terrain) on each model level. 
#
# Dependencies  : netCDF4, numpy, scipy
#
  
  
from __future__ import print_function # make sure print behaves the same in Python 2.7 and 3.x
import netCDF4
from netCDF4 import num2date, Dataset
import numpy as np
from scipy import interpolate
import datetime
import sys
import io
import math
import pandas as pd
import os

#Routine to Read File A data file for t, q values at a particular grid point
def readfa(fileapath):
    #Connect to data file for reading
    fa = netCDF4.Dataset(fileapath, 'r')
  
    fanclongs = fa.variables['longitude']  
    fanclats = fa.variables['latitude']  
    fanclevels = fa.variables['level']  
    fanctimes = fa.variables['time']  
    fancts = fa.variables['t'] 
    fancqs = fa.variables['q'] 
  
    #Get level values (either model levels or pressure levels) and number of levels
    falevels=fanclevels[:]
    fanlevels=falevels.shape[0]
  
    #Get time values and number of times
    fatimes=fanctimes[:]
    fantimes=fatimes.shape[0]
  
    #Get python datetime for each time
    fapydts=num2date(fatimes, fanctimes.units)
  
    #Get lat and long values
    falats = fanclats[:]
    falongs = fanclongs[:]
  
    #Get index of grid point of interest
    failat=np.where(falats==GRID_LAT)[0]
    failong=np.where(falongs==GRID_LONG)[0]
  
    #Get t, q values for specified grid point for all levels (slicing)
    #The result is still a 4D array with 1 latitude and 1 longitude
    fats=fancts[range(fantimes),range(fanlevels),failat,failong]
    faqs=fancqs[range(fantimes),range(fanlevels),failat,failong]
  
    return (fats, faqs, falevels, fapydts)
  
#Routine to Read File B data file for z, lnsp values at a particular grid point
def readfb(fbfilepath):
    #Connect to file for reading
    fbf = netCDF4.Dataset(fbfilepath, 'r')
  
    fbnclongs = fbf.variables['longitude']  
    fbnclats = fbf.variables['latitude']  
    fbnctimes = fbf.variables['time']  
    fbnczs = fbf.variables['z'] 
    fbnclnsps = fbf.variables['lnsp'] 
  
    #Get time values and number of times
    fbtimes=fbnctimes[:]
    fbntimes=fbtimes.shape[0]
  
    #Get python datetime for each time
    fbpydts=num2date(fbtimes, fbnctimes.units)
  
    #Output first 10 datetimes
  
    #Get lat and long values
    fblats = fbnclats[:]
    fblongs = fbnclongs[:]
  
  
    #Get index of grid point of interest
    fbilat=np.where(fblats==GRID_LAT)[0]
    fbilong=np.where(fblongs==GRID_LONG)[0]
  
    #Get z, lnsp values for specified grid point (slicing)
    #The result is still a 3D array with 1 latitude and 1 longitude
    fbzs=fbnczs[range(fbntimes),fbilat,fbilong]
    fblnsps=fbnclnsps[range(fbntimes),fbilat,fbilong]
  
    return (fbzs, fblnsps, fbpydts)

def calculatepress(z, lnsp, ts, qs, levels):
    presreturn=np.full(ts.shape[0], -999, np.double)
      
    #surface pressure
    sp = math.exp(lnsp)
  
    # A and B parameters to calculate pressures for model levels, 
    #  extracted from an ECMWF ERA-Interim GRIB file and then hardcoded here
    pv =  [0,2.000365,3.102241,4.666084,6.827977,9.746966,13.605424,18.608931,24.985718,32.98571,42.879242,54.955463,69.520576,86.895882,107.415741,131.425507,
    159.279404,191.338562,227.968948,269.539581,316.420746,368.982361,427.592499,492.616028,564.413452,643.339905,729.744141,823.967834,926.34491,1037.201172,
    1156.853638,1285.610352,1423.770142,1571.622925,1729.448975,1897.519287,2076.095947,2265.431641,2465.770508,2677.348145,2900.391357,3135.119385,3381.743652,
    3640.468262,3911.490479,4194.930664,4490.817383,4799.149414,5119.89502,5452.990723,5798.344727,6156.074219,6526.946777,6911.870605,7311.869141,7727.412109,
    8159.354004,8608.525391,9076.400391,9562.682617,10065.97852,10584.63184,11116.66211,11660.06738,12211.54785,12766.87305,13324.66895,13881.33106,14432.13965,
    14975.61523,15508.25684,16026.11523,16527.32227,17008.78906,17467.61328,17901.62109,18308.43359,18685.71875,19031.28906,19343.51172,19620.04297,19859.39063,
    20059.93164,20219.66406,20337.86328,20412.30859,20442.07813,20425.71875,20361.81641,20249.51172,20087.08594,19874.02539,19608.57227,19290.22656,18917.46094,
    18489.70703,18006.92578,17471.83984,16888.6875,16262.04688,15596.69531,14898.45313,14173.32422,13427.76953,12668.25781,11901.33984,11133.30469,10370.17578,
    9617.515625,8880.453125,8163.375,7470.34375,6804.421875,6168.53125,5564.382813,4993.796875,4457.375,3955.960938,3489.234375,3057.265625,2659.140625,
    2294.242188,1961.5,1659.476563,1387.546875,1143.25,926.507813,734.992188,568.0625,424.414063,302.476563,202.484375,122.101563,62.78125,22.835938,3.757813,
    0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0.000007,0.000024,0.000059,0.000112,0.000199,0.00034,0.000562,0.00089,0.001353,0.001992,0.002857,0.003971,0.005378,0.007133,0.009261,0.011806,
    0.014816,0.018318,0.022355,0.026964,0.032176,0.038026,0.044548,0.051773,0.059728,0.068448,0.077958,0.088286,0.099462,0.111505,0.124448,0.138313,
    0.153125,0.16891,0.185689,0.203491,0.222333,0.242244,0.263242,0.285354,0.308598,0.332939,0.358254,0.384363,0.411125,0.438391,0.466003,0.4938,
    0.521619,0.549301,0.576692,0.603648,0.630036,0.655736,0.680643,0.704669,0.727739,0.749797,0.770798,0.790717,0.809536,0.827256,0.843881,0.859432,
    0.873929,0.887408,0.8999,0.911448,0.922096,0.931881,0.94086,0.949064,0.95655,0.963352,0.969513,0.975078,0.980072,0.984542,0.9885,0.991984,0.995003,
    0.99763,1
      ]
    levelSize=137
    A = pv[0:levelSize+1]
    B = pv[levelSize+1:]
  
    Ph_levplusone = A[levelSize] + (B[levelSize]*sp)
  
    #Get a list of level numbers in reverse order
    reversedlevels=np.full(levels.shape[0], -999, np.int32)
    for iLev in list(reversed(range(levels.shape[0]))):
        reversedlevels[levels.shape[0] - 1 - iLev] = levels[iLev]
              
    #Integrate up into the atmosphere from lowest level
    for lev in reversedlevels:  
        #compute the pressures (on half-levels)
        ilevel=np.where(levels==lev)[0]
        
        if lev == 137:
            presreturn[ilevel] = sp
            
            Ph_lev = sp
            
        else:

            Ph_lev = A[lev-1] + (B[lev-1] * sp)
            
            presreturn[ilevel] = (Ph_levplusone + Ph_lev)/2
   
        Ph_levplusone = Ph_lev
  
    return presreturn

df = pd.read_csv(r'C:\Users\dylan\Downloads\NHP Hail Work\Hail_data_MATEUSZ\Chosen Profiles.csv')
df['datetime'] = pd.to_datetime(df['ERA5_DATE'])

for k in df['datetime']:
    file_name = k.strftime('%Y-%m-%d-%H')
    if os.path.isfile(r'C:\Users\dylan\Downloads\NHP Hail Work\ERA5 Complete Data/' + f'ERA5-{file_name}-p.nc'):
        pass
    else:

        FILE_A_PATH=r"C:\Users\dylan\Downloads\NHP Hail Work\ERA5 Complete Data/"+f"ERA5-{file_name}-tq_ml.grib" 
        FILE_B_PATH=r"C:\Users\dylan\Downloads\NHP Hail Work\ERA5 Complete Data/"+f"ERA5-{file_name}-zlnp.grib"
        
        fa = netCDF4.Dataset(FILE_A_PATH, 'r')
        fb = netCDF4.Dataset(FILE_B_PATH, 'r')
        
        ncfile = Dataset(r'C:\Users\dylan\Downloads\NHP Hail Work\ERA5 Complete Data/'+f'ERA5-{file_name}-p.nc',mode='w',format='NETCDF4_CLASSIC') 
        
        lat_dim = ncfile.createDimension('lat', len(fb['latitude'][:]))
        lon_dim = ncfile.createDimension('lon', len(fb['longitude'][:]))
        level_dim = ncfile.createDimension('level', len(fa['level'][:]))
        time_dim = ncfile.createDimension('time', 1)
        ncfile.title='Model Level Pressures'
        
        lat = ncfile.createVariable('lat', np.float32, ('lat',))
        lat.units = 'degrees_north'
        lat.long_name = 'latitude'
        lon = ncfile.createVariable('lon', np.float32, ('lon',))
        lon.units = 'degrees_east'
        lon.long_name = 'longitude'
        time = ncfile.createVariable('time', np.float64, ('time',))
        time.units = 'hours since 1800-01-01'
        time.long_name = 'time'
        level = ncfile.createVariable('level', np.float32, ('level',))
        
        pres = ncfile.createVariable('p',np.float64,('time','level','lat','lon')) # note: unlimited dimension is leftmost
        pres.units = 'Pa' # degrees Kelvin
        pres.standard_name = 'pressure' # this is a CF standard name
        
        lat[:] = fb['latitude'][:]
        lon[:] = fb['longitude'][:]
        
        for i in range(len(fb['latitude'][:])):
            for j in range(len(fb['longitude'][:])):
                GRID_LAT=fb['latitude'][i]   #(degrees N)
                GRID_LONG=fb['longitude'][j] #(degrees W)     
                  
                #Read File A file 
                fats, faqs, falevels, fapydts = readfa(FILE_A_PATH)
                fanlevels=falevels.shape[0]
                fantimes=fapydts.shape[0]
                  
                #Read File B file 
                fbzs, fblnsps, fbpydts = readfb(FILE_B_PATH)
                      
                unicode = str
                
                #Output all the values for the specified grid point
                #Iterate over all the times and write out the data
                for itime in range(fantimes):
                                
                    pydt=fbpydts[itime]
                    z, lnsp = fbzs[itime,0,0], fblnsps[itime,0,0]
                    sline="{}, {}, {}".format(pydt.strftime('%Y-%m-%d %H:%M:%S'), z, lnsp)
                    
                    ps = calculatepress(z, lnsp, fats[itime,range(fanlevels),0,0],
                        faqs[itime,range(fanlevels),0,0], falevels)
                    
                    pres[:,:,i,j] = ps
        
        ncfile.close()       
