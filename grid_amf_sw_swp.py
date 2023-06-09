#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May  9 22:36:54 2023

@author: benjaminyang
"""
##########################################################################################
#####   Grid Aura-OMI L2 air mass factors, scattering weights, and pressure levels   #####
##########################################################################################

# NOTE: This script may take ~65 hours to run (depends on system).  
# Recommended to run in the background (e.g. use "screen" command in Linux). 

# Import modules
import numpy as np
from datetime import datetime, date, time, timedelta
from dateutil.relativedelta import relativedelta
from netCDF4 import Dataset
import glob

# Get all dates between start and end dates inclusive
def date_range(start, end, delta):
    curr = start
    while curr <= end:
        yield curr
        curr += delta

# Read latitude and longitude from file into numpy arrays
def naive_fast(latvar,lonvar,lat0,lon0):
    latvals = latvar[:]
    lonvals = lonvar[:]
    ny,nx = latvals.shape
    dist_sq = (latvals-lat0)**2 + (lonvals-lon0)**2
    minindex_flattened = dist_sq.argmin()  # 1D index of min element
    iy_min,ix_min = np.unravel_index(minindex_flattened, latvals.shape)
    return int(iy_min),int(ix_min)    

# Enter start and end dates
start = date(2018,1,1)
end = date(2018,12,31)

# Save all-inclusive dates to a list
alldates=[]
for result in date_range(start, end, relativedelta(months=1)):
    alldates.append(result.strftime('%Y%m'))
alldates

# Define paths of data and figure directories
fig_dir = '/home/byang/output/amore/' # figures
obs_dir = '/home/byang/data/' # observational data

# Get 1D arrays of lat/lon from OMI L3 NO2 and HCHO files
nc0 = Dataset(obs_dir+'OMI_VCD/EasternUS_OMI_monthlymean_NO2_2018.nc','r')
lat_no2 = nc0.variables['lat'][:] # 0.25 deg
lon_no2 = nc0.variables['lon'][:] # 0.25 deg
nc0.close()
nc0 = Dataset(obs_dir+'OMI_VCD/EasternUS_OMI_monthlymean_HCHO_2018.nc','r')
lat_hcho = nc0.variables['lat'][:] # 0.1 deg
lon_hcho = nc0.variables['lon'][:] # 0.1 deg
nc0.close()

# Save all-inclusive dates to a list
alldates2=[]
for result in date_range(start, end, relativedelta(days=1)):
    alldates2.append(result.strftime('%m%d'))
alldates2

# Build 2D grids (rows = lat, cols = lon) for AMF_m, SW_m, and SWP_m
amf_grid_all = [] # month=12,species=2,lat=108/270,lon=144/360
sw_grid_all = [] # month=12,species=2,lat=108/270,lon=144/360
swp_grid_all = [] # month=12,species=1,lat=270,lon=360 (HCHO only)
species = ['NO2','HCHO']
slat = [lat_no2,lat_hcho]
slon = [lon_no2,lon_hcho]
slev = [35,47] # number of vertical pressure levels [NO2,HCHO]
for d in alldates: # loop over months in 2018
    #print(d)
    amf_grid_m = [[],[]]
    sw_grid_m = [[],[]]
    swp_grid_m = []
    mdates = [d2 for d2 in alldates2 if d2.startswith(d[4:])]   
    for s,ss in enumerate(species): # loop over each species
        #dd = mdates[0]
        for dd in mdates: # loop over days in given month
            #print(dd)
            rows = len(slat[s])
            columns= len(slon[s])
            amf_grid = [[np.nan for x in range(columns)] for x in range(rows)]
            sw_grid = [[[np.nan]*slev[s] for x in range(columns)] for x in range(rows)]
            swp_grid = [[[np.nan]*slev[s] for x in range(columns)] for x in range(rows)]

            # Create sorted list of OMI swath files for given day
            omi_path = sorted(glob.glob('/data0/byang/OMI_VCD/OMI_%s_L2/OMI-Aura_L2-OM%s_2018m%s*'%(ss,ss,dd)))
            if not omi_path: # if files are missing for given day (skip)
                print('Missing '+dd)
                continue

            lat = []
            lon = []
            amf = []
            sw = []
            swp = []
            for o in omi_path: # loop over swaths for given day
                # Read in data variables
                hdf = Dataset(o,'r')
                if ss=='NO2':
                    var = hdf.groups['HDFEOS'].groups['SWATHS'].groups['ColumnAmountNO2']
                    lat.append(var.groups['Geolocation Fields']['Latitude'][:]) # Latitude of the center of the groundpixel (degrees N) (dimensions: nTimes=1644, nXtrack=60)
                    lon.append(var.groups['Geolocation Fields']['Longitude'][:]) # Longitude of the center of the groundpixel (degrees E) (dimensions: nTimes=1644, nXtrack=60)
                    amf.append(var.groups['Data Fields']['AmfTrop'][:]) # Tropospheric air mass factor (dimensions: nTimes=1644, nXtrack=60)
                    sw.append(var.groups['Data Fields']['ScatteringWeight'][:]) # Scattering weight profile (nTimes=1644, nXtrack=60, nPresLevels=35)
                elif ss=='HCHO':
                    var = hdf.groups['HDFEOS'].groups['SWATHS'].groups['OMI Total Column Amount HCHO']
                    lat.append(var.groups['Geolocation Fields']['Latitude'][:]) # Geodetic Latitude (degrees N) (dimensions: nTimes=1644, nXtrack=60)
                    lon.append(var.groups['Geolocation Fields']['Longitude'][:]) # Geodetic Longitude (degrees E) (dimensions: nTimes=1644, nXtrack=60)
                    amf.append(var.groups['Data Fields']['AirMassFactor'][:]) # Molecule Specific Air Mass Factor (AMF) (dimensions: nTimes=1643, nXtrack=60)
                    sw.append(np.transpose(var.groups['Data Fields']['ScatteringWeights'][:],(1,2,0))) # Scattering Weights (nLevels=47, nTimes=1643, nXtrack=60)
                    swp.append(np.transpose(var.groups['Data Fields']['ClimatologyLevels'][:],(1,2,0))) # Climatology Levels (hPa) (nLevels=47, nTimes=1643, nXtrack=60)
                hdf.close()

            # Reshape 3d arrays into 2d (combine all swaths for given day)
            try: # if all swaths have same shape (faster)
                lat = np.reshape(lat,(1644*len(omi_path),60))
                lon = np.reshape(lon,(1644*len(omi_path),60))
                amf = np.reshape(amf,(1644*len(omi_path),60))
                sw = np.reshape(sw,(1644*len(omi_path),60,slev[s]))
            except: # if swaths have different shapes (slower)
                lat = np.concatenate(lat,axis=0)
                lon = np.concatenate(lon,axis=0)
                amf = np.concatenate(amf,axis=0)
                sw = np.concatenate(sw,axis=0)
            if ss=='HCHO':
                try:
                    swp = np.reshape(swp,(1644*len(omi_path),60,slev[s]))
                except:
                    swp = np.concatenate(swp,axis=0)

            for i in range(rows):
                for j in range(columns):    
                    # Find the nearest lat/lon point (of all swaths) to the given grid box
                    iy,ix = naive_fast(lat,lon,slat[s][i],slon[s][j])

                    # Add AMF_m and SW_m to grids
                    if amf[iy,ix]>=0: # valid value
                        amf_grid[i][j] = amf[iy,ix]
                    else: # if value is masked (large negative value)
                        continue
                    if sw[iy,ix,0]>=0: # valid value
                        sw_grid[i][j] = sw[iy,ix]
                    else: # if value is masked (large negative value)
                        continue
                    if ss=='HCHO':
                        if swp[iy,ix,0]>=0: # valid value
                            swp_grid[i][j] = swp[iy,ix]
                        else: # if value is masked (large negative value)
                            continue

                    #print(j)
                #print(i)
                print('date: 2018%s, var: %s, lat: %d, lon: %d'%(dd,ss,i,j))

            amf_grid_m[s].append(amf_grid)
            sw_grid_m[s].append(sw_grid)
            if ss=='HCHO':
                swp_grid_m.append(swp_grid)

    amf_grid_all.append([np.nanmean(amf_grid_m[i],axis=0) for i in [0,1]])
    sw_grid_all.append([np.nanmean(sw_grid_m[i],axis=0) for i in [0,1]])
    swp_grid_all.append(np.nanmean(swp_grid_m,axis=0))

# Save converted data arrays to numpy files
np.save(fig_dir+'amf_grid_all.npy',np.array(amf_grid_all))
np.save(fig_dir+'sw_grid_all.npy',np.array(sw_grid_all))
np.save(fig_dir+'swp_grid_all.npy',np.array(swp_grid_all))
