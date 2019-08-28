#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 23 14:30:56 2019

@author: ricardofaria

Script to create Large-scale forcing data netcdf input for PALM following PALM Input Data Standard (PIDS) v1.10.
my_test_setup_radiation.nc – contains all radiation data need to force PALM radiation.
Contains static and dynamic information of radiation properties (trace gas profiles, sky view factors).

"""


case_name = 'madeira_50m_offline_nesting'   # PALM case name to fit WRF boundary conditions data
wrf_domain = 'd02_airport'
zone = '28'
interp_mode = 'linear'  # linear
#wrf_geos_dom_buffer = 15 # Grid buffer for WRF domain boundary's geostrophic profiles calc, number of cells from boundarys


import numpy as np
from netCDF4 import Dataset, num2date
from wrf import getvar, ALL_TIMES 
import glob
import pandas as pd
import res_grid_change as rgc
import nearest


cfg = pd.read_csv(case_name + '.cfg')
dx = cfg.dx.values[0]
dy = cfg.dy.values[0]
dz = cfg.dz.values[0]
nx = cfg.nx.values[0]
ny = cfg.ny.values[0]
nz = cfg.nz.values[0]
north = cfg.north.values[0]
south = cfg.south.values[0]
east = cfg.east.values[0]
west = cfg.west.values[0]

nc_st = Dataset(case_name + '_static', 'r')
y, x, z = nc_st.variables['y'][:], nc_st.variables['x'][:], np.arange(dz/2, dz*nz, dz)
xu = x + np.gradient(x)/2
xu = xu[:-1]
yv = y + np.gradient(y)/2
yv = yv[:-1]
zw = z + np.gradient(z)/2
zw = zw[:-1]
zt = nc_st.variables['zt'][:]
zt = dz * np.round(zt/dz)
origin_lat = nc_st.origin_lat
origin_lon = nc_st.origin_lon
nc_st.close()


file_forcing = sorted(glob.glob('raw_forcing/*' + wrf_domain + '*')) 
times = []
for ff_idx, ff in enumerate(file_forcing) :
    print('Loading WRF netCDF: ' + str(ff))
    nc_ff = Dataset(ff, 'r')
    if ff_idx == 0 :
        lat_wrf = nc_ff.variables['XLAT'][0,:,0].data
        lon_wrf = nc_ff.variables['XLONG'][0,0,:].data
        
        south_idx, north_idx = nearest.nearest(lat_wrf, south)[1], nearest.nearest(lat_wrf, north)[1]
        west_idx, east_idx = nearest.nearest(lon_wrf, west)[1], nearest.nearest(lon_wrf, east)[1]
        
        lat_wrf = nc_ff.variables['XLAT'][south_idx:north_idx].data
        lon_wrf = nc_ff.variables['XLONG'][west_idx:east_idx].data
        
        z_wrf = getvar(nc_ff, 'height', units='m')[:,south_idx:north_idx, west_idx:east_idx] # nc_ff.variables['height'][:]
        zstag_wrf = getvar(nc_ff, 'zstag', units='m')[:,south_idx:north_idx, west_idx:east_idx] # nc_ff.variables['zstag'][:]
        
        zs = nc_ff.variables['ZS'][0, :].data
        dzs = nc_ff.variables['DZS'][0, :].data
        landmask = nc_ff.variables['LANDMASK'][0, south_idx:north_idx, west_idx:east_idx].data
        tmn = nc_ff.variables['TMN'][0, south_idx:north_idx, west_idx:east_idx].data
        
    
    if ff_idx == 0 :
        rad_sw_in = getvar(nc_ff, 'SWDOWN', timeidx = ALL_TIMES)[:, south_idx:north_idx, west_idx:east_idx]
        rad_lw_in = getvar(nc_ff, 'GLW', timeidx = ALL_TIMES)[:, south_idx:north_idx, west_idx:east_idx]
    else :
        rad_sw_in = np.append(rad_sw_in, np.atleast_3d(getvar(nc_ff, 'SWDOWN', timeidx = ALL_TIMES)[:, south_idx:north_idx, west_idx:east_idx]), axis=0)
        rad_lw_in = np.append(rad_lw_in, np.atleast_3d(getvar(nc_ff, 'GLW', timeidx = ALL_TIMES)[:, south_idx:north_idx, west_idx:east_idx]), axis=0)
    
    time_var = nc_ff.variables['XTIME']
    times = np.append(times, num2date(time_var[:],time_var.units))
    
    nc_ff.close()

tot_sec = (times[-1]-times[0]).total_seconds()
time_setp_sec = (times[1]-times[0]).total_seconds()
times_sec = np.arange(0, tot_sec+1, time_setp_sec)


rad_sw_in_tmp = np.zeros((rad_sw_in.shape[0], y.shape[0], x.shape[0]))
rad_lw_in_tmp = np.zeros((rad_lw_in.shape[0], y.shape[0], x.shape[0]))


for t in range(rad_sw_in_tmp.shape[0]):
    
    rad_sw_in_tmp[t, :, :] = rgc.grid_points_change(rad_sw_in[t, :, :], x.shape[0], y.shape[0], interp_mode)
    rad_lw_in_tmp[t, :, :] = rgc.grid_points_change(rad_lw_in[t, :, :], x.shape[0], y.shape[0], interp_mode)
    
    

nc_output = Dataset(case_name + '_dynamic', 'r+')
# Global Attributes
#nc_output.description = 'Contains radiation data from WRF mesoscale'
#nc_output.history = 'Created by Ricardo Faria - OOM (ricardo88faria@gmail.com) ' + time.ctime(time.time())
#nc_output.source = 'netCDF4 python'
#nc_output.source = 'netCDF4 python'
#nc_output.origin_lat = np.float(origin_lat)
#nc_output.origin_lon = np.float(origin_lon)
#nc_output.origin_z = np.float(0)
#nc_output.origin_x = np.float(0)
#nc_output.origin_y = np.float(0)
#nc_output.rotation_angle = np.float(0)
#nc_output.origin_time = str(times[0]) + ' +00'

#nc_output.createDimension('x', nx)
#nc_output.createDimension('y', ny)
nc_output.createDimension('time_rad', len(times_sec))

#nc_x = nc_output.createVariable('x',  np.float32, 'x')
#nc_y = nc_output.createVariable('y',  np.float32, 'y')
nc_time_rad = nc_output.createVariable('time_rad', np.float32, 'time_rad')

nc_rad_sw_in = nc_output.createVariable('rad_sw_in',  np.float32, ('time_rad', 'y', 'x'), fill_value = -9999.0, zlib=True)
nc_rad_lw_in = nc_output.createVariable('rad_lw_in',  np.float32, ('time_rad', 'y', 'x'), fill_value = -9999.0, zlib=True)

#nc_x.units = 'm'
#nc_y.units = 'm'
nc_time_rad.units = 'seconds'

nc_rad_sw_in.long_name = 'incoming shortwave radiative flux at the surfacev (W/m^2)'
nc_rad_sw_in.source = 'WRF'
nc_rad_sw_in.units = 'W/m^2'
nc_rad_sw_in.lod = np.float(2)
nc_rad_lw_in.long_name = 'incoming longwave radiative flux at the surface (W/m^2)'
nc_rad_lw_in.source = 'WRF'
nc_rad_lw_in.units = 'W/m^2'
nc_rad_lw_in.lod = np.float(2)

#nc_x[:] = x
#nc_y[:] = y
nc_time_rad[:] = times_sec

rad_sw_in_tmp[rad_sw_in_tmp == np.nan] = 0
rad_lw_in_tmp[rad_lw_in_tmp == np.nan] = 0

nc_rad_sw_in[:] = rad_sw_in_tmp
nc_rad_lw_in[:] = rad_lw_in_tmp

nc_output.close()

