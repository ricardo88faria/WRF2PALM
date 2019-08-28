#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul  2 12:08:53 2018

@author: ricardofaria

Script to create Large-scale forcing data (boundary conditions) netcdf input for PALM following PALM Input Data Standard (PIDS) v1.9.
my_test_setup_dynamic.nc – contains all dynamic data need to force PALM.
Contains dynamic information for the run, such as time­dependent boundary conditions and the initial state of the atmosphere.

"""


case_name = 'madeira_breeze_50m_offline_nesting'   # PALM case name to fit WRF boundary conditions data
wrf_domain = 'd02_breeze'
zone = '28'
interp_mode = 'linear'  # linear
wrf_geos_dom_buffer = 15 # Grid buffer for WRF domain boundary's geostrophic profiles calc, number of cells from boundarys


import numpy as np
from netCDF4 import Dataset, num2date
from wrf import getvar, ALL_TIMES, interplevel#, vinterp, destagger
import time
import glob
import pandas as pd
from pyproj import Proj 
import scipy.interpolate as interpolate
import res_grid_change as rgc
import nearest
import geostrophic


dz_soil = np.array([0.01, 0.02, 0.04, 0.06, 0.14, 0.26, 0.54, 1.86])

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
nzt_ns = np.array([zt[0, :], zt[-1, :]])/dz
nzt_ew = np.array([zt[:, 0], zt[:, -1]])/dz
origin_lat = nc_st.origin_lat
origin_lon = nc_st.origin_lon
nc_st.close()

z_int_lev = np.arange(z[0], z[-1] + nz*dz + 1, dz)

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
        
    
    U = getvar(nc_ff, 'ua', timeidx = ALL_TIMES)[:,:,south_idx:north_idx, west_idx:east_idx] # nc_ff.variables['U'][:].data
    V = getvar(nc_ff, 'va', timeidx = ALL_TIMES)[:,:,south_idx:north_idx, west_idx:east_idx] # nc_ff.variables['V'][:].data
    W = getvar(nc_ff, 'wa', timeidx = ALL_TIMES)[:,:,south_idx:north_idx, west_idx:east_idx] # nc_ff.variables['W'][:].data
    qvapor = nc_ff.variables['QVAPOR'][:,:,south_idx:north_idx, west_idx:east_idx].data
    theta = getvar(nc_ff, 'theta', timeidx = ALL_TIMES, units='K')[:,:,south_idx:north_idx, west_idx:east_idx]
    TSLB = nc_ff.variables['TSLB'][:,:,south_idx:north_idx, west_idx:east_idx].data
    SMOIS = nc_ff.variables['SMOIS'][:,:,south_idx:north_idx, west_idx:east_idx].data
    
    
    qvapor_int = np.empty((qvapor.shape[0], len(z_int_lev), qvapor.shape[2], qvapor.shape[3]))
    theta_int = np.empty((qvapor.shape[0], len(z_int_lev), qvapor.shape[2], qvapor.shape[3]))
    U_int = np.empty((qvapor.shape[0], len(z_int_lev), qvapor.shape[2], qvapor.shape[3]))
    V_int = np.empty((qvapor.shape[0], len(z_int_lev), qvapor.shape[2], qvapor.shape[3]))
    W_int = np.empty((qvapor.shape[0], len(z_int_lev), qvapor.shape[2], qvapor.shape[3]))
    
    # specific by z horizontal levels interpolation
    for l_idx, l in enumerate(z_int_lev) :
        print('Interpolating height level: ' + str(l) + ' [m]')
        if l == 0 :
            l = 1
            z_int = zstag_wrf[:-1]
        else :
            z_int = z_wrf
        
#        qvapor_interp = interplevel(qvapor, z_int, l)
        qvapor_int[:, int(l_idx), :] = interplevel(qvapor, z_int, l).data
        
#        theta_interp = interplevel(theta, z_int, l)
        theta_int[:, int(l_idx), :] = interplevel(theta, z_int, l).data
        
#        U_interp = interplevel(U, z_int, l)
        U_int[:, int(l_idx), :] = interplevel(U, z_int, l).data
        
#        V_interp = interplevel(V, z_int, l)
        V_int[:, int(l_idx), :] = interplevel(V, z_int, l).data
        
#        W_interp = interplevel(W, z_wrf, l)
        W_int[:, int(l_idx), :] = interplevel(W, z_wrf, l).data
    
    if ff_idx == 0 :
        u = U_int #u_tot = U_int
        v = V_int
        w = W_int
        qv = qvapor_int
        pt = theta_int

        tslb = TSLB
        smois = SMOIS
        
        
    else :
        u = np.append(u, np.atleast_3d(U_int), axis=0)
        v = np.append(v, np.atleast_3d(V_int), axis=0)
        w = np.append(w, np.atleast_3d(W_int), axis=0)
        qv = np.append(qv, np.atleast_3d(qvapor_int), axis=0)
        pt = np.append(pt, np.atleast_3d(theta_int), axis=0)
        tslb = np.append(tslb, np.atleast_3d(TSLB), axis=0)
        smois = np.append(smois, np.atleast_3d(SMOIS), axis=0)

    time_var = nc_ff.variables['XTIME']
    times = np.append(times, num2date(time_var[:],time_var.units))
    
    del time_var, U_int, V_int, W_int, qvapor_int, theta_int, U, V, W, qvapor, theta, TSLB, SMOIS,  #, T2, Q2, PSFC
    
    nc_ff.close()

tot_sec = (times[-1]-times[0]).total_seconds()
time_setp_sec = (times[1]-times[0]).total_seconds()
times_sec = np.arange(0, tot_sec+1, time_setp_sec)

#change latlon projection to UTM meters
#myProj = Proj("+proj=utm +zone=" + zone + ", +north +ellps=WGS84 +datum=WGS84 +units=m +no_defs")
#inProj = Proj(init='epsg:4326')
#outProj = Proj(init='epsg:3061')

# create north-south, east-west and top boundary conditions and interpolate to new xyz shape
u_ns_bc_tmp = np.zeros((times.shape[0], z.shape[0], 2, x.shape[0]))
u_ew_bc_tmp = np.zeros((times.shape[0], z.shape[0], y.shape[0], 2))
u_top_bc_tmp = np.zeros((times.shape[0], y.shape[0], x.shape[0]))

v_ns_bc_tmp = np.zeros((times.shape[0], z.shape[0], 2, x.shape[0]))
v_ew_bc_tmp = np.zeros((times.shape[0], z.shape[0], y.shape[0], 2))
v_top_bc_tmp = np.zeros((times.shape[0], y.shape[0], x.shape[0]))

w_ns_bc_tmp = np.zeros((times.shape[0], z.shape[0], 2, x.shape[0]))
w_ew_bc_tmp = np.zeros((times.shape[0], z.shape[0], y.shape[0], 2))
w_top_bc_tmp = np.zeros((times.shape[0], y.shape[0], x.shape[0]))

qv_ns_bc_tmp = np.zeros((times.shape[0], z.shape[0], 2, x.shape[0]))
qv_ew_bc_tmp = np.zeros((times.shape[0], z.shape[0], y .shape[0], 2))
qv_top_bc_tmp = np.zeros((times.shape[0], y.shape[0], x.shape[0]))

pt_ns_bc_tmp = np.zeros((times.shape[0], z.shape[0], 2, x.shape[0]))
pt_ew_bc_tmp = np.zeros((times.shape[0], z.shape[0], y.shape[0], 2))
pt_top_bc_tmp = np.zeros((times.shape[0], y.shape[0], x.shape[0]))

idx = [0, -1]
for t in range(u.shape[0]):
    
#     bring surface to 0 index to set height above ground for PALM as  1st step for better interpolation values
    for x_idx in range(u.shape[3]) :
        for y_idx in range(u.shape[2]) :
            
            nan_idx = np.argwhere(np.isnan(u[t,:,y_idx,x_idx]))
            
            if nan_idx.size != 0 :
                nan_idx = nan_idx[-1] + 1
                
                u[t,:,y_idx,x_idx] = rgc.grid_points_change_1darray(u[t,int(nan_idx):,y_idx,x_idx], u.shape[1])
                v[t,:,y_idx,x_idx] = rgc.grid_points_change_1darray(v[t,int(nan_idx):,y_idx,x_idx], v.shape[1])
                w[t,:,y_idx,x_idx] = rgc.grid_points_change_1darray(w[t,int(nan_idx):,y_idx,x_idx], w.shape[1])
                pt[t,:,y_idx,x_idx] = rgc.grid_points_change_1darray(pt[t,int(nan_idx):,y_idx,x_idx], pt.shape[1])
                qv[t,:,y_idx,x_idx] = rgc.grid_points_change_1darray(qv[t,int(nan_idx):,y_idx,x_idx], qv.shape[1])
    
    for i in idx:
        
        
        u_ns_bc_tmp[t, :, i , :] = rgc.grid_points_change(u[t, :, i, :], x.shape[0], z.shape[0], interp_mode)
        v_ns_bc_tmp[t, :, i , :] = rgc.grid_points_change(v[t, :, i, :], x.shape[0], z.shape[0], interp_mode)
        w_ns_bc_tmp[t, :, i , :] = rgc.grid_points_change(w[t, :, i, :], x.shape[0], z.shape[0], interp_mode)
        qv_ns_bc_tmp[t, :, i , :] = rgc.grid_points_change(qv[t, :, i, :], x.shape[0], z.shape[0], interp_mode)
        pt_ns_bc_tmp[t, :, i , :] = rgc.grid_points_change(pt[t, :, i, :], x.shape[0], z.shape[0], interp_mode)
        
        for x_idx in range(u_ns_bc_tmp.shape[3]) :
            
            u_ns_bc_tmp[t, int(nzt_ns[i, x_idx]):, i, x_idx] = rgc.grid_points_change_1darray(u_ns_bc_tmp[t, :, i, x_idx], nz - nzt_ns[i, x_idx])
            u_ns_bc_tmp[t, :int(nzt_ns[i, x_idx]), i, x_idx] = np.nan
            v_ns_bc_tmp[t, int(nzt_ns[i, x_idx]):, i, x_idx] = rgc.grid_points_change_1darray(v_ns_bc_tmp[t, :, i, x_idx], nz - nzt_ns[i, x_idx])
            v_ns_bc_tmp[t, :int(nzt_ns[i, x_idx]), i, x_idx] = np.nan
            w_ns_bc_tmp[t, int(nzt_ns[i, x_idx]):, i, x_idx] = rgc.grid_points_change_1darray(w_ns_bc_tmp[t, :, i, x_idx], nz - nzt_ns[i, x_idx])
            w_ns_bc_tmp[t, :int(nzt_ns[i, x_idx]), i, x_idx] = np.nan
            qv_ns_bc_tmp[t, int(nzt_ns[i, x_idx]):, i, x_idx] = rgc.grid_points_change_1darray(qv_ns_bc_tmp[t, :, i, x_idx], nz - nzt_ns[i, x_idx])
            qv_ns_bc_tmp[t, :int(nzt_ns[i, x_idx]), i, x_idx] = np.nan
            pt_ns_bc_tmp[t, int(nzt_ns[i, x_idx]):, i, x_idx] = rgc.grid_points_change_1darray(pt_ns_bc_tmp[t, :, i, x_idx], nz - nzt_ns[i, x_idx])
            pt_ns_bc_tmp[t, :int(nzt_ns[i, x_idx]), i, x_idx] = np.nan
        
        u_ew_bc_tmp[t, :, : , i] = rgc.grid_points_change(u[t, :, :, i], y.shape[0], z.shape[0], interp_mode)
        v_ew_bc_tmp[t, :, : , i] = rgc.grid_points_change(v[t, :, :, i], y.shape[0], z.shape[0], interp_mode)
        w_ew_bc_tmp[t, :, : , i] = rgc.grid_points_change(w[t, :, :, i], y.shape[0], z.shape[0], interp_mode)
        qv_ew_bc_tmp[t, :, : , i] = rgc.grid_points_change(qv[t, :, :, i], y.shape[0], z.shape[0], interp_mode)
        pt_ew_bc_tmp[t, :, : , i] = rgc.grid_points_change(pt[t, :, :, i], y.shape[0], z.shape[0], interp_mode)
        
        for y_idx in range(u_ew_bc_tmp.shape[2]) :
                
            u_ew_bc_tmp[t, int(nzt_ew[i, y_idx]):, y_idx, i] = rgc.grid_points_change_1darray(u_ew_bc_tmp[t, :, y_idx, i], nz - nzt_ew[i, y_idx])
            u_ew_bc_tmp[t, :int(nzt_ew[i, y_idx]), y_idx, i] = np.nan
            v_ew_bc_tmp[t, int(nzt_ew[i, y_idx]):, y_idx, i] = rgc.grid_points_change_1darray(v_ew_bc_tmp[t, :, y_idx, i], nz - nzt_ew[i, y_idx])
            v_ew_bc_tmp[t, :int(nzt_ew[i, y_idx]), y_idx, i] = np.nan
            w_ew_bc_tmp[t, int(nzt_ew[i, y_idx]):, y_idx, i] = rgc.grid_points_change_1darray(w_ew_bc_tmp[t, :, y_idx, i], nz - nzt_ew[i, y_idx])
            w_ew_bc_tmp[t, :int(nzt_ew[i, y_idx]), y_idx, i] = np.nan
            qv_ew_bc_tmp[t, int(nzt_ew[i, y_idx]):, y_idx, i] = rgc.grid_points_change_1darray(qv_ew_bc_tmp[t, :, y_idx, i], nz - nzt_ew[i, y_idx])
            qv_ew_bc_tmp[t, :int(nzt_ew[i, y_idx]), y_idx, i] = np.nan
            pt_ew_bc_tmp[t, int(nzt_ew[i, y_idx]):, y_idx, i] = rgc.grid_points_change_1darray(pt_ew_bc_tmp[t, :, y_idx, i], nz - nzt_ew[i, y_idx])
            pt_ew_bc_tmp[t, :int(nzt_ew[i, y_idx]), y_idx, i] = np.nan
#        
    u_top_bc_tmp[t, :, :] = rgc.grid_points_change(u[t, -1, :, :], x.shape[0], y.shape[0], interp_mode)
    v_top_bc_tmp[t, :, :] = rgc.grid_points_change(v[t, -1, :, :], x.shape[0], y.shape[0], interp_mode)
    w_top_bc_tmp[t, :, :] = rgc.grid_points_change(w[t, -1, :, :], x.shape[0], y.shape[0], interp_mode)
    qv_top_bc_tmp[t, :, :] = rgc.grid_points_change(qv[t, -1, :, :], x.shape[0], y.shape[0], interp_mode)
    pt_top_bc_tmp[t, :, :] = rgc.grid_points_change(pt[t, -1, :, :], x.shape[0], y.shape[0], interp_mode)

## geostrophic profiles calc
file_forcing_d01 = sorted(glob.glob('raw_forcing/*' + wrf_domain + '*')) 
for ff_idx, ff in enumerate(file_forcing_d01) :
    nc_ff = Dataset(ff, 'r')
    
    if ff_idx == 0 :
        lat_wrf = nc_ff.variables['XLAT'][0,:,0].data
        lon_wrf = nc_ff.variables['XLONG'][0,0,:].data
        z_wrf = getvar(nc_ff, 'height', units='m')[:]
#        zstag_wrf = getvar(nc_ff, 'zstag', units='m')[:]
        
    tk = getvar(nc_ff, 'tk', timeidx = ALL_TIMES)
    phid = getvar(nc_ff, 'pres', timeidx = ALL_TIMES, units='Pa')
    
    if ff_idx == 0 :
        
        temp_tot = tk
        ph_tot = phid
        
    else :
        temp_tot = np.append(temp_tot, np.atleast_3d(tk), axis=0)
        ph_tot = np.append(ph_tot, np.atleast_3d(phid), axis=0)


lat_wrf_f, lon_wrf_f = lat_wrf[wrf_geos_dom_buffer:-wrf_geos_dom_buffer], lon_wrf[wrf_geos_dom_buffer:-wrf_geos_dom_buffer]

# still to implement the point (lat lon) where your profile should be calculated !!!
geo_wind_u = np.empty((ph_tot.shape[0], ph_tot.shape[1]))
geo_wind_v = np.empty((ph_tot.shape[0], ph_tot.shape[1]))
geo_wind_u_f = np.empty((u.shape[0], u.shape[1]))
geo_wind_v_f = np.empty((u.shape[0], u.shape[1]))
for t in range(ph_tot.shape[0]) :
    for h in range(1, ph_tot.shape[1]) :
        geo_wind = geostrophic.geostr(ph_tot[t, h, wrf_geos_dom_buffer:-wrf_geos_dom_buffer], temp_tot[t, h, wrf_geos_dom_buffer:-wrf_geos_dom_buffer], lat_wrf_f[:], lon_wrf_f[:])
        geo_wind_u[t, h] = geo_wind[0]
        geo_wind_v[t, h] = geo_wind[1]
    
    geo_wind_u[t, :][np.logical_and( geo_wind_u[t, :] >= -.1, geo_wind_u[t, :] <= .1)] = np.nan
    geo_wind_u[t, :][np.logical_or( geo_wind_u[t, :] >= 99, geo_wind_u[t, :] <= -99)] = np.nan
    geo_wind_v[t, :][np.logical_and( geo_wind_v[t, :] >= -.1, geo_wind_v[t, :] <= .1)] = np.nan
    geo_wind_v[t, :][np.logical_or( geo_wind_v[t, :] >= 99, geo_wind_v[t, :] <= -99)] = np.nan
    
    A = geo_wind_u[t, :]
    inds = np.arange(A.shape[0])
    good = np.where(np.isfinite(A))
    f = interpolate.interp1d(inds[good], A[good], fill_value='extrapolate')
    geo_wind_u[t, :] = np.where(np.isfinite(A), A, f(inds))
    A = geo_wind_v[t, :]
    inds = np.arange(A.shape[0])
    good = np.where(np.isfinite(A))
    f = interpolate.interp1d(inds[good], A[good], fill_value='extrapolate')
    geo_wind_v[t, :] = np.where(np.isfinite(A), A, f(inds))
    
#    plt.plot(geo_wind_u[t, :], z_wrf.data[:, 0, 0])
#    plt.plot(geo_wind_v[t, :], z_wrf.data[:, 0, 0])
    
    geo_wind_u_f[t, :] = np.interp(z, z_wrf.data[:, 0, 0], geo_wind_u[t, :])
    geo_wind_v_f[t, :] = np.interp(z, z_wrf.data[:, 0, 0], geo_wind_v[t, :])
    
#    plt.plot(geo_wind_u_f[t, :], z)
#    plt.plot(geo_wind_v_f[t, :], z)
#    plt.show()
    


nc_output = Dataset(case_name + '_dynamic', 'w')
# Global Attributes
nc_output.description = 'Contains dynamic data from WRF mesoscale'
nc_output.history = 'Created by Ricardo Faria - OOM (ricardo88faria@gmail.com) ' + time.ctime(time.time())
nc_output.source = 'netCDF4 python'
nc_output.source = 'netCDF4 python'
nc_output.origin_lat = np.float(-17)
nc_output.origin_lon = np.float(32)
nc_output.origin_z = np.float(0)
nc_output.origin_x = np.float(0)
nc_output.origin_y = np.float(0)
nc_output.rotation_angle = np.float(0)
nc_output.origin_time = str(times[0]) + ' +00'

nc_output.createDimension('x', nx)
nc_output.createDimension('y', ny)
nc_output.createDimension('z', nz)
nc_output.createDimension('zsoil', len(dz_soil))
nc_output.createDimension('xu', nx-1)
nc_output.createDimension('yv', ny-1)
nc_output.createDimension('zw', nz-1)
nc_output.createDimension('time', len(times_sec))

nc_x = nc_output.createVariable('x',  np.float32, 'x')
nc_y = nc_output.createVariable('y',  np.float32, 'y')  
nc_z = nc_output.createVariable('z', np.float32, 'z')
nc_zsoil = nc_output.createVariable('zsoil', np.float32, 'zsoil')
nc_xu = nc_output.createVariable('xu',  np.float32, 'xu')
nc_yv = nc_output.createVariable('yv',  np.float32, 'yv')  
nc_zw = nc_output.createVariable('zw', np.float32, 'zw')
nc_time = nc_output.createVariable('time', np.float32, 'time')

nc_init_soil_m = nc_output.createVariable('init_soil_m',  np.float32, ('zsoil'), fill_value = -9999.0, zlib=True)
nc_init_soil_t = nc_output.createVariable('init_soil_t',  np.float32, ('zsoil'), fill_value = -9999.0, zlib=True)

nc_ls_forcing_ug = nc_output.createVariable('ls_forcing_ug',  np.float32, ('time', 'z'), fill_value = -9999.0, zlib=True)
nc_ls_forcing_vg = nc_output.createVariable('ls_forcing_vg',  np.float32, ('time', 'z'), fill_value = -9999.0, zlib=True)

nc_ls_forcing_left_u = nc_output.createVariable('ls_forcing_left_u',  np.float32, ('time', 'z', 'y'), fill_value = -9999.0, zlib=True)
nc_ls_forcing_left_v = nc_output.createVariable('ls_forcing_left_v',  np.float32, ('time', 'z', 'yv'), fill_value = -9999.0, zlib=True)
nc_ls_forcing_left_w = nc_output.createVariable('ls_forcing_left_w',  np.float32, ('time', 'zw', 'y'), fill_value = -9999.0, zlib=True)
nc_ls_forcing_left_qv = nc_output.createVariable('ls_forcing_left_qv',  np.float32, ('time', 'z', 'y'), fill_value = -9999.0, zlib=True)
nc_ls_forcing_left_pt = nc_output.createVariable('ls_forcing_left_pt',  np.float32, ('time', 'z', 'y'), fill_value = -9999.0, zlib=True)

nc_ls_forcing_right_u = nc_output.createVariable('ls_forcing_right_u',  np.float32, ('time', 'z', 'y'), fill_value = -9999.0, zlib=True)
nc_ls_forcing_right_v = nc_output.createVariable('ls_forcing_right_v',  np.float32, ('time', 'z', 'yv'), fill_value = -9999.0, zlib=True)
nc_ls_forcing_right_w = nc_output.createVariable('ls_forcing_right_w',  np.float32, ('time', 'zw', 'y'), fill_value = -9999.0, zlib=True)
nc_ls_forcing_right_qv = nc_output.createVariable('ls_forcing_right_qv',  np.float32, ('time', 'z', 'y'), fill_value = -9999.0, zlib=True)
nc_ls_forcing_right_pt = nc_output.createVariable('ls_forcing_right_pt',  np.float32, ('time', 'z', 'y'), fill_value = -9999.0, zlib=True)

nc_ls_forcing_north_u = nc_output.createVariable('ls_forcing_north_u',  np.float32, ('time', 'z', 'xu'), fill_value = -9999.0, zlib=True)
nc_ls_forcing_north_v = nc_output.createVariable('ls_forcing_north_v',  np.float32, ('time', 'z', 'x'), fill_value = -9999.0, zlib=True)
nc_ls_forcing_north_w = nc_output.createVariable('ls_forcing_north_w',  np.float32, ('time', 'zw', 'x'), fill_value = -9999.0, zlib=True)
nc_ls_forcing_north_qv = nc_output.createVariable('ls_forcing_north_qv',  np.float32, ('time', 'z', 'x'), fill_value = -9999.0, zlib=True)
nc_ls_forcing_north_pt = nc_output.createVariable('ls_forcing_north_pt',  np.float32, ('time', 'z', 'x'), fill_value = -9999.0, zlib=True)

nc_ls_forcing_south_u = nc_output.createVariable('ls_forcing_south_u',  np.float32, ('time', 'z', 'xu'), fill_value = -9999.0, zlib=True)
nc_ls_forcing_south_v = nc_output.createVariable('ls_forcing_south_v',  np.float32, ('time', 'z', 'x'), fill_value = -9999.0, zlib=True)
nc_ls_forcing_south_w = nc_output.createVariable('ls_forcing_south_w',  np.float32, ('time', 'zw', 'x'), fill_value = -9999.0, zlib=True)
nc_ls_forcing_south_qv = nc_output.createVariable('ls_forcing_south_qv',  np.float32, ('time', 'z', 'x'), fill_value = -9999.0, zlib=True)
nc_ls_forcing_south_pt = nc_output.createVariable('ls_forcing_south_pt',  np.float32, ('time', 'z', 'x'), fill_value = -9999.0, zlib=True)

nc_ls_forcing_top_u = nc_output.createVariable('ls_forcing_top_u',  np.float32, ('time', 'y', 'xu'), fill_value = -9999.0, zlib=True)
nc_ls_forcing_top_v = nc_output.createVariable('ls_forcing_top_v',  np.float32, ('time', 'yv', 'x'), fill_value = -9999.0, zlib=True)
nc_ls_forcing_top_w = nc_output.createVariable('ls_forcing_top_w',  np.float32, ('time', 'y', 'x'), fill_value = -9999.0, zlib=True)
nc_ls_forcing_top_qv = nc_output.createVariable('ls_forcing_top_qv',  np.float32, ('time', 'y', 'x'), fill_value = -9999.0, zlib=True)
nc_ls_forcing_top_pt = nc_output.createVariable('ls_forcing_top_pt',  np.float32, ('time', 'y', 'x'), fill_value = -9999.0, zlib=True)


nc_x.units = 'm'
nc_y.units = 'm'
nc_z.units = 'm'
nc_zsoil.units = 'm'
nc_xu.units = 'm'
nc_yv.units = 'm'
nc_zw.units = 'm'
nc_time.units = 'seconds'

nc_init_soil_m.long_name = 'volumetric soil moisture (m^3/m^3)'
nc_init_soil_m.source = 'WRF'
nc_init_soil_m.units = 'm^3/m^3'
nc_init_soil_m.lod = np.float(1)
nc_init_soil_t.long_name = 'soil temperature (K)'
nc_init_soil_t.source = 'WRF'
nc_init_soil_t.units = 'K'
nc_init_soil_t.lod = np.float(1)

nc_ls_forcing_ug.long_name = 'u wind component geostrophic'
nc_ls_forcing_ug.source = 'WRF'
nc_ls_forcing_ug.units = 'm/s'
nc_ls_forcing_ug.standard_name = ''
nc_ls_forcing_vg.long_name = 'v wind component geostrophic'
nc_ls_forcing_vg.source = 'WRF'
nc_ls_forcing_vg.units = 'm/s'
nc_ls_forcing_vg.standard_name = ''

nc_ls_forcing_left_u.long_name = 'ls_forcing_left_u'
nc_ls_forcing_left_u.res_orig = str(dx) + 'x' + str(dy) + ' m'
nc_ls_forcing_left_u.source = 'WRF'
nc_ls_forcing_left_u.units = 'm/s'
nc_ls_forcing_left_u.standard_name = ''
nc_ls_forcing_left_v.long_name = 'ls_forcing_left_v'
nc_ls_forcing_left_v.res_orig = str(dx) + 'x' + str(dy) + ' m'
nc_ls_forcing_left_v.source = 'WRF'
nc_ls_forcing_left_v.units = 'm/s'
nc_ls_forcing_left_v.standard_name = ''
nc_ls_forcing_left_w.long_name = 'ls_forcing_left_w'
nc_ls_forcing_left_w.res_orig = str(dx) + 'x' + str(dy) + ' m'
nc_ls_forcing_left_w.source = 'WRF'
nc_ls_forcing_left_w.units = 'm/s'
nc_ls_forcing_left_w.standard_name = ''
nc_ls_forcing_left_qv.long_name = 'ls_forcing_left_qv'
nc_ls_forcing_left_qv.res_orig = str(dx) + 'x' + str(dy) + ' m'
nc_ls_forcing_left_qv.source = 'WRF'
nc_ls_forcing_left_qv.units = 'kg/kg'
nc_ls_forcing_left_qv.standard_name = ''
nc_ls_forcing_left_pt.long_name = 'ls_forcing_left_pt'
nc_ls_forcing_left_pt.res_orig = str(dx) + 'x' + str(dy) + ' m'
nc_ls_forcing_left_pt.source = 'WRF'
nc_ls_forcing_left_pt.units = 'K'
nc_ls_forcing_left_pt.standard_name = ''

nc_ls_forcing_right_u.long_name = 'ls_forcing_right_u'
nc_ls_forcing_right_u.res_orig = str(dx) + 'x' + str(dy) + ' m'
nc_ls_forcing_right_u.source = 'WRF'
nc_ls_forcing_right_u.units = 'm/s'
nc_ls_forcing_right_u.standard_name = ''
nc_ls_forcing_right_v.long_name = 'ls_forcing_right_v'
nc_ls_forcing_right_v.res_orig = str(dx) + 'x' + str(dy) + ' m'
nc_ls_forcing_right_v.source = 'WRF'
nc_ls_forcing_right_v.units = 'm/s'
nc_ls_forcing_right_v.standard_name = ''
nc_ls_forcing_right_w.long_name = 'ls_forcing_right_w'
nc_ls_forcing_right_w.res_orig = str(dx) + 'x' + str(dy) + ' m'
nc_ls_forcing_right_w.source = 'WRF'
nc_ls_forcing_right_w.units = 'm/s'
nc_ls_forcing_right_w.standard_name = ''
nc_ls_forcing_right_qv.long_name = 'ls_forcing_right_qv'
nc_ls_forcing_right_qv.res_orig = str(dx) + 'x' + str(dy) + ' m'
nc_ls_forcing_right_qv.source = 'WRF'
nc_ls_forcing_right_qv.units = 'kg/kg'
nc_ls_forcing_right_qv.lod = np.float(2)
nc_ls_forcing_right_pt.long_name = 'ls_forcing_right_pt'
nc_ls_forcing_right_pt.res_orig = str(dx) + 'x' + str(dy) + ' m'
nc_ls_forcing_right_pt.source = 'WRF'
nc_ls_forcing_right_pt.units = 'K'
nc_ls_forcing_right_pt.standard_name = ''

nc_ls_forcing_north_u.long_name = 'ls_forcing_north_u'
nc_ls_forcing_north_u.res_orig = str(dx) + 'x' + str(dy) + ' m'
nc_ls_forcing_north_u.source = 'WRF'
nc_ls_forcing_north_u.units = 'm/s'
nc_ls_forcing_north_u.standard_name = ''
nc_ls_forcing_north_v.long_name = 'ls_forcing_north_v'
nc_ls_forcing_north_v.res_orig = str(dx) + 'x' + str(dy) + ' m'
nc_ls_forcing_north_v.source = 'WRF'
nc_ls_forcing_north_v.units = 'm/s'
nc_ls_forcing_north_v.standard_name = ''
nc_ls_forcing_north_w.long_name = 'ls_forcing_north_w'
nc_ls_forcing_north_w.res_orig = str(dx) + 'x' + str(dy) + ' m'
nc_ls_forcing_north_w.source = 'WRF'
nc_ls_forcing_north_w.units = 'm/s'
nc_ls_forcing_north_w.standard_name = ''
nc_ls_forcing_north_qv.long_name = 'ls_forcing_north_qv'
nc_ls_forcing_north_qv.res_orig = str(dx) + 'x' + str(dy) + ' m'
nc_ls_forcing_north_qv.source = 'WRF'
nc_ls_forcing_north_qv.units = 'kg/kg'
nc_ls_forcing_north_qv.lod = np.float(2)
nc_ls_forcing_north_qv.standard_name = ''
nc_ls_forcing_north_pt.long_name = 'ls_forcing_north_pt'
nc_ls_forcing_north_pt.res_orig = str(dx) + 'x' + str(dy) + ' m'
nc_ls_forcing_north_pt.source = 'WRF'
nc_ls_forcing_north_pt.units = 'K'
nc_ls_forcing_north_pt.standard_name = ''

nc_ls_forcing_south_u.long_name = 'ls_forcing_south_u'
nc_ls_forcing_south_u.res_orig = str(dx) + 'x' + str(dy) + ' m'
nc_ls_forcing_south_u.source = 'WRF'
nc_ls_forcing_south_u.units = 'm/s'
nc_ls_forcing_south_u.standard_name = ''
nc_ls_forcing_south_v.long_name = 'ls_forcing_south_v'
nc_ls_forcing_south_v.res_orig = str(dx) + 'x' + str(dy) + ' m'
nc_ls_forcing_south_v.source = 'WRF'
nc_ls_forcing_south_v.units = 'm/s'
nc_ls_forcing_south_v.standard_name = ''
nc_ls_forcing_south_w.long_name = 'ls_forcing_south_w'
nc_ls_forcing_south_w.res_orig = str(dx) + 'x' + str(dy) + ' m'
nc_ls_forcing_south_w.source = 'WRF'
nc_ls_forcing_south_w.units = 'm/s'
nc_ls_forcing_south_w.standard_name = ''
nc_ls_forcing_south_qv.long_name = 'ls_forcing_south_qv'
nc_ls_forcing_south_qv.res_orig = str(dx) + 'x' + str(dy) + ' m'
nc_ls_forcing_south_qv.source = 'WRF'
nc_ls_forcing_south_qv.units = 'kg/kg'
nc_ls_forcing_south_qv.standard_name = ''
nc_ls_forcing_south_pt.long_name = 'ls_forcing_south_pt'
nc_ls_forcing_south_pt.res_orig = str(dx) + 'x' + str(dy) + ' m'
nc_ls_forcing_south_pt.source = 'WRF'
nc_ls_forcing_south_pt.units = 'K'
nc_ls_forcing_south_pt.standard_name = ''

nc_ls_forcing_top_u.long_name = 'ls_forcing_top_u'
nc_ls_forcing_top_u.res_orig = str(dx) + 'x' + str(dy) + ' m'
nc_ls_forcing_top_u.source = 'WRF'
nc_ls_forcing_top_u.units = 'm/s'
nc_ls_forcing_top_u.standard_name = ''
nc_ls_forcing_top_v.long_name = 'ls_forcing_top_v'
nc_ls_forcing_top_v.res_orig = str(dx) + 'x' + str(dy) + ' m'
nc_ls_forcing_top_v.source = 'WRF'
nc_ls_forcing_top_v.units = 'm/s'
nc_ls_forcing_top_w.long_name = 'ls_forcing_top_w'
nc_ls_forcing_top_w.res_orig = str(dx) + 'x' + str(dy) + ' m'
nc_ls_forcing_top_w.source = 'WRF'
nc_ls_forcing_top_w.units = 'm/s'
nc_ls_forcing_top_qv.long_name = 'ls_forcing_top_qv'
nc_ls_forcing_top_qv.res_orig = str(dx) + 'x' + str(dy) + ' m'
nc_ls_forcing_top_qv.source = 'WRF'
nc_ls_forcing_top_qv.units = 'kg/kg'
nc_ls_forcing_top_pt.long_name = 'ls_forcing_top_pt'
nc_ls_forcing_top_pt.res_orig = str(dx) + 'x' + str(dy) + ' m'
nc_ls_forcing_top_pt.source = 'WRF'
nc_ls_forcing_top_pt.units = 'K'


nc_x[:] = x #np.arange(0, dx*nx, dx)
nc_y[:] = y #np.arange(0, dy*ny, dy)
nc_z[:] = z
nc_zsoil[:] = dz_soil
nc_xu[:] = xu
nc_yv[:] = yv
nc_zw[:] = zw
nc_time[:] = times_sec

nc_ls_forcing_ug[:] = geo_wind_u_f
nc_ls_forcing_vg[:] = geo_wind_v_f

top_u = np.empty(nc_ls_forcing_top_u.shape)
top_v = np.empty(nc_ls_forcing_top_v.shape)
top_w = np.empty(nc_ls_forcing_top_w.shape)
top_qv = np.empty(nc_ls_forcing_top_qv.shape)
top_pt = np.empty(nc_ls_forcing_top_pt.shape)


init_soil_t = np.empty(zs.shape)
init_soil_m = np.empty(zs.shape)
init_soil_tmn = np.nanmean(np.ma.masked_where(landmask == 0, tmn))

for d in range(zs.shape[0]) :
    
    init_soil_t[d] = np.nanmean(np.ma.masked_where(landmask == 0, tslb[0,d,:,:]))
    init_soil_m[d] = np.nanmean(np.ma.masked_where(landmask == 0, smois[0,d,:,:]))


init_soil_t = np.interp(dz_soil, zs, init_soil_t)
init_soil_m = np.interp(dz_soil, zs, init_soil_m)

nc_init_soil_m[:] = init_soil_m
nc_init_soil_t[:] = init_soil_t

u_ns_bc_tmp[u_ns_bc_tmp == np.nan] = 0
u_ew_bc_tmp[u_ew_bc_tmp == np.nan] = 0
v_ns_bc_tmp[v_ns_bc_tmp == np.nan] = 0
v_ew_bc_tmp[v_ew_bc_tmp == np.nan] = 0
w_ns_bc_tmp[w_ns_bc_tmp == np.nan] = 0
w_ew_bc_tmp[w_ew_bc_tmp == np.nan] = 0
qv_ns_bc_tmp[qv_ns_bc_tmp == np.nan] = 0
qv_ew_bc_tmp[qv_ew_bc_tmp == np.nan] = 0
pt_ns_bc_tmp[pt_ns_bc_tmp == np.nan] = 0
pt_ew_bc_tmp[pt_ew_bc_tmp == np.nan] = 0

nc_ls_forcing_left_u[:] = u_ew_bc_tmp[:, :, :, 0]
nc_ls_forcing_left_v[:] = v_ew_bc_tmp[:, :, :-1, 0]
nc_ls_forcing_left_w[:] = w_ew_bc_tmp[:, :-1, :, 0]
nc_ls_forcing_left_qv[:] = qv_ew_bc_tmp[:, :, :, 0]
nc_ls_forcing_left_pt[:] = pt_ew_bc_tmp[:, :, :, 0]

nc_ls_forcing_right_u[:] = u_ew_bc_tmp[:, :, :, 1]
nc_ls_forcing_right_v[:] = v_ew_bc_tmp[:, :, :-1, 1]
nc_ls_forcing_right_w[:] = w_ew_bc_tmp[:, :-1, :, 1]
nc_ls_forcing_right_qv[:] = qv_ew_bc_tmp[:, :, :, 1]
nc_ls_forcing_right_pt[:] = pt_ew_bc_tmp[:, :, :, 1]

nc_ls_forcing_north_u[:] = u_ns_bc_tmp[:, :, 1, :-1]
nc_ls_forcing_north_v[:] = v_ns_bc_tmp[:, :, 1, :]
nc_ls_forcing_north_w[:] = w_ns_bc_tmp[:, :-1, 1, :]
nc_ls_forcing_north_qv[:] = qv_ns_bc_tmp[:, :, 1, :]
nc_ls_forcing_north_pt[:] = pt_ns_bc_tmp[:, :, 1, :]

nc_ls_forcing_south_u[:] = u_ns_bc_tmp[:, :, 0, :-1]
nc_ls_forcing_south_v[:] = v_ns_bc_tmp[:, :, 0, :]
nc_ls_forcing_south_w[:] = w_ns_bc_tmp[:, :-1, 0, :]
nc_ls_forcing_south_qv[:] = qv_ns_bc_tmp[:, :, 0, :]
nc_ls_forcing_south_pt[:] = pt_ns_bc_tmp[:, :, 0, :]

nc_ls_forcing_top_u[:] = u_top_bc_tmp[:,:,:-1]
nc_ls_forcing_top_v[:] = v_top_bc_tmp[:,:-1,:]
nc_ls_forcing_top_w[:] = w_top_bc_tmp
nc_ls_forcing_top_qv[:] = qv_top_bc_tmp
nc_ls_forcing_top_pt[:] = pt_top_bc_tmp

nc_output.close()

print('Add to your *_p3d file the: ' + '\n soil_temperature = ' + repr(init_soil_t) + '\n soil_moisture = ' + repr(init_soil_m) + '\n deep_soil_temperature = ' + repr(init_soil_tmn))

