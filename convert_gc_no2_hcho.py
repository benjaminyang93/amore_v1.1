#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jan 21 22:49:57 2023

@author: benjaminyang
"""

#############################################################################################################################
#####   Convert GEOS-Chem NO2 and HCHO mole fractions to vertical column densities + calculate model air mass factors   #####
#############################################################################################################################

# NOTE: This script may take ~42 hours to run (depends on system).
# Recommended to run in the background (e.g. use "screen" command in Linux). 
# Must have already run "grid_amf_sw_swp.py" first.

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

# Enter start and end dates
start = date(2018,1,1)
end = date(2018,12,31)

# Save all-inclusive dates to a list
alldates=[]
for result in date_range(start, end, relativedelta(months=1)):
    alldates.append(result.strftime('%Y%m'))
alldates

# Define paths of data and figure directories
data_dir = '/data0/byang/amore_v4/' # AMORE
data_dir2 = '/data0/byang/gc_merra2_fullchem_complexSOA_SVPOA/' # BASE
data_dir3 = '/data0/byang/gc_merra2_fullchem_complexSOA_SVPOA_zUS_ISOP/' # BASE zUS_ISOP 
data_dir4 = '/data0/byang/amore_zUS_ISOP/' # AMORE zUS_ISOP 
data_dir5 = '/data0/byang/gc_merra2_fullchem_complexSOA_SVPOA_zUS_ANTH/' # BASE zUS_ANTH 
data_dir6 = '/data0/byang/amore_zUS_ANTH/' # AMORE zUS_ANTH 
data_dir7 = '/data0/byang/gc_merra2_fullchem_complexSOA_SVPOA_zUS_ISOP_ANTH/' # BASE zUS_ISOP_ANTH 
data_dir8 = '/data0/byang/amore_zUS_ISOP_ANTH/' # AMORE zUS_ISOP_ANTH 
data_dir9 = '/data0/byang/gc_merra2_fullchem_complexSOA_SVPOA_zGLB_ANTH/' # BASE zGLB_ANTH 
data_dir10 = '/data0/byang/amore_zGLB_ANTH/' # AMORE zGLB_ANTH 
data_dir11 = '/data0/byang/gc_geosfp_fullchem_complexSOA_SVPOA_US/' # BASE_HR
data_dir12 = '/data0/byang/amore_geosfp_highres/' # AMORE_HR
fig_dir = '/home/byang/output/amore/' # figures
obs_dir = '/home/byang/data/' # observational data
all_data_dir = [data_dir,data_dir2,data_dir3,data_dir4,data_dir5,data_dir6,data_dir7,data_dir8,data_dir9,data_dir10,data_dir11,data_dir12] # all 12 GEOS-Chem simulations

# Get 1D arrays of lat/lon from OMI L3 NO2 and HCHO files
nc0 = Dataset(obs_dir+'OMI_VCD/EasternUS_OMI_monthlymean_NO2_2018.nc','r')
lat_no2 = nc0.variables['lat'][:] # 0.25 deg
lon_no2 = nc0.variables['lon'][:] # 0.25 deg
nc0.close()
nc0 = Dataset(obs_dir+'OMI_VCD/EasternUS_OMI_monthlymean_HCHO_2018.nc','r')
lat_hcho = nc0.variables['lat'][:] # 0.1 deg
lon_hcho = nc0.variables['lon'][:] # 0.1 deg
nc0.close()

# Read in saved data arrays from numpy files
sw_grid_all = np.load(fig_dir+'sw_grid_all.npy',allow_pickle=True) # month=12,species=2,lat=108/270,lon=144/360,lev=35/47
swp_grid_all = np.load(fig_dir+'swp_grid_all.npy',allow_pickle=True) # month=12,lat=270,lon=360,lev=47 (HCHO only)

# Read in fixed scattering weight pressures from an NO2 HDF data file
hdf = Dataset('/data0/byang/OMI_VCD/OMI_NO2_L2/OMI-Aura_L2-OMNO2_2018m0101t1530-o71624_v003-2019m0819t193019.he5','r')
var = hdf.groups['HDFEOS'].groups['SWATHS'].groups['ColumnAmountNO2']
swp = var.groups['Data Fields']['ScatteringWtPressure'][:] # Pressures for scattering weight profile (hPa) (nPresLevels=35)
hdf.close()

# Convert model NO2 and HCHO using MERRA2 (2° x 2.5°) and GEOSFP (0.25° x 0.3125°) data
vcd = [] # month=12,species=2,simulation=10,lat=91,lon=144 (coarse)
vcd_hr = [] # month=12,species=2,simulation=2,lat=100,lon=108 (high-res)
amf = [] # month=12,species=2,simulation=10,lat=91,lon=144 (coarse)
amf_hr = [] # month=12,species=2,simulation=2,lat=100,lon=108 (high-res)
for di,d in enumerate(alldates):
    print('Starting...')
    for s in [10,2]: # coarser (10) and higher resolution (2) sensitivity simulations  
        if s==10:
            ps = 'C' # coarse-res
            met_path = sorted(glob.glob('/data0/byang/MERRA2/MERRA2.%s*.A1.2x25.nc4'%d))
            met_path2 = sorted(glob.glob('/data0/byang/MERRA2/MERRA2.%s*.I3.2x25.nc4'%d))
        else:
            ps = 'H' # high-res
            met_path = sorted(glob.glob('/data0/byang/GEOSFP/GEOSFP.%s*.A1.025x03125.nc'%d))
            met_path2 = sorted(glob.glob('/data0/byang/GEOSFP/GEOSFP.%s*.I3.025x03125.nc'%d))

        # Read in MERRA2 or GEOSFP meteorological data and average daily then monthly    
        tropp = []
        psfc = []
        T = []
        q = []
        for m,m2 in zip(met_path,met_path2):
            met = Dataset(m,'r')
            met2 = Dataset(m2,'r')
            if s==10:
                tropp.append(np.mean(met.variables['TROPPT'],axis=0)) # tropopause pressure based on thermal estimate (Pa)
                psfc.append(np.mean(met2.variables['PS'],axis=0)) # surface pessure (Pa)
            else:
                tropp.append(np.mean(met.variables['TROPPT'],axis=0)*100) # tropopause pressure based on thermal estimate (hPa -> Pa)
                psfc.append(np.mean(met2.variables['PS'],axis=0)*100) # surface pessure (hPa -> Pa)
            T.append(np.mean(met2.variables['T'],axis=0)) # air temperature (K)
            q.append(np.mean(met2.variables['QV'],axis=0)) # specific humidity (kg/kg)
            met.close()
            met2.close()
            #print(m)
        tropp_m = np.mean(tropp,axis=0)
        psfc_m = np.mean(psfc,axis=0)
        T_m = np.mean(T,axis=0)
        q_m = np.mean(q,axis=0)

        # Read in and combine GEOS-Chem model output for all simulations (coarser or higher-res)
        sn = 'GEOSChem.SpeciesConc.%s01_0000z.nc4'%d
        nc = []
        for n in range(s):
            if s==2:
                n = n+10
            nc.append(Dataset(all_data_dir[n]+sn,'r')) 

        # Convert model NO2 and HCHO (ppb) to tropospheric VCD (molec/cm^2) for each grid box monthly
        species = ['SpeciesConc_NO2','SpeciesConc_CH2O']
        shape = np.shape(nc[0].variables[species[0]][:]) 
        var = [[[[0 for ilon in range(shape[3])] for ilat in range(shape[2])] for n in range(s)] for v in species]
        var_new = [[[[0 for ilon in range(shape[3])] for ilat in range(shape[2])] for n in range(s)] for v in species]
        for ilat in range(shape[2]): # loop over latitudes 
            for ilon in range(shape[3]): # loop over longitudes 
                var_vcd = [[[] for n in range(s)] for v in species]
                var_vcd_new = [[[] for n in range(s)] for v in species]
                for l in range(shape[1]): # loop over vertical levels 
                    ilev1 = nc[0].variables['ilev'][l] # hybrid level at interfaces ("bottom edge")
                    ilev2 = nc[0].variables['ilev'][l+1] # hybrid level at interfaces ("top edge")
                    lev = nc[0].variables['lev'][l] # hybrid level at midpoints ((A/P0)+B)
                    lat = nc[0].variables['lat'][ilat] # Latitude (degrees N)
                    lon = nc[0].variables['lon'][ilon] # Longitude (degrees E)

                    # Get meteorological variables for grid box
                    T = T_m[l,ilat,ilon] # air temperature (K)
                    q = q_m[l,ilat,ilon] # specific humidity (kg/kg)

                    # Compute pressures for grid box
                    P = lev*psfc_m[ilat,ilon] # pressure at center of grid box (Pa), where P_top_atmosphere = 0
                    P1 = ilev1*psfc_m[ilat,ilon] # pressure at bottom of grid box (Pa), where P_top_atmosphere = 0
                    P2 = ilev2*psfc_m[ilat,ilon] # pressure at top of grid box (Pa), where P_top_atmosphere = 0

                    # Define constants
                    R = 8.3145 # universal gas constant (J/(K*mol), J = Pa*m^3)
                    AN = 6.02214*1e23 # Avogadro's number (molecules/mol)
                    g = 9.81 # gravitational acceleration (m/s^2)
                    Rd = 287 # specific gas constant for dry air (J/(kg*K))

                    # Solve hysometric equation
                    w = q/(1-q) # water vapor mixing ratio (kg/kg)
                    Tv = T*(1+(0.608*w)) # aprrox. virtual temperature (K)
                    h = Rd*Tv/g*np.log(P1/P2) # thickness of layer (m)

                    # Species NO2 or HCHO (ppb)
                    for iv,v in enumerate(species):     
                        for n in range(s):
                            var_ppb = nc[n].variables[v][0,l,ilat,ilon]*1e9 

                            # Convert model NO2 or HCHO from ppb to molec/cm^2 
                            var_molec_cm2 = var_ppb/1e9*P/(R*T)*AN*h/1e4 # use ideal gas law (n/V = p/RT)
                            
                            if v=='SpeciesConc_CH2O': # HCHO: entire vertical column
                                var_vcd[iv][n].append(var_molec_cm2)
                                # Multiply partial column by scattering weight
                                ix = np.argmin(np.abs(lat_hcho - lat))
                                iy = np.argmin(np.abs(lon_hcho - lon))
                                iz = np.argmin(np.abs(swp_grid_all[di,ix,iy,:] - P/100))
                                var_vcd_new[iv][n].append(var_molec_cm2*sw_grid_all[di,1][ix,iy,iz])
                            else: # NO2: tropospheric vertical column
                                if P2 > tropp_m[ilat,ilon]: # below the tropopause 
                                    var_vcd[iv][n].append(var_molec_cm2)
                                    # Multiply partial column by scattering weight
                                    ix = np.argmin(np.abs(lat_no2 - lat))
                                    iy = np.argmin(np.abs(lon_no2 - lon))
                                    iz = np.argmin(np.abs(swp - P/100))
                                    var_vcd_new[iv][n].append(var_molec_cm2*sw_grid_all[di,0][ix,iy,iz])
                                else: # at or above the tropopause 
                                    break

                # Integrate (add up) levels in each tropospheric (NO2) or entire (HCHO) vertical column 
                for iv,v in enumerate(species):
                    for n in range(s):
                        var[iv][n][ilat][ilon] = sum(var_vcd[iv][n]) # VCD_m
                        var_new[iv][n][ilat][ilon] = sum(var_vcd_new[iv][n])/var[iv][n][ilat][ilon] # AMF_m
                print('date: %s, res: %s, lat: %d, lon: %d'%(d,ps,ilat,ilon))

        if s==10:
            vcd.append(var)
            amf.append(var_new)
        else:
            vcd_hr.append(var)
            amf_hr.append(var_new)

# Save converted data arrays to numpy files
np.save(fig_dir+'vcd.npy',np.array(vcd))
np.save(fig_dir+'vcd_hr.npy',np.array(vcd_hr))
np.save(fig_dir+'amf.npy',np.array(amf))
np.save(fig_dir+'amf_hr.npy',np.array(amf_hr))

