#!/usr/bin/env python3
# -*- coding: utf-8 -*-§
"""

Created on Tue Jun 26 14:13:28 2018

@author: ricardofaria

Script to create static driver netcdf input for PALM following PALM Input Data Standard (PIDS) v1.9
my_test_setup_static.nc – contains all static information like orography, buildings, and surface classification.

Still to implement corine land cover information.

"""


# user input
topo_file_name = 'raw_topo/srtm_33_06.asc'
case_name = 'madeira_breeze_50m_offline_nesting'
zone = '28'
north = 32.89
south = 32.61
east = -16.64
west = -17.29
dx = 50
dy = 50
dz = 50
nz = 80



import numpy as np
from netCDF4 import Dataset
import time
import pandas as pd
from pyproj import Proj
import read_geo as rg
import res_grid_change as rgc
import nearest
import matplotlib.pyplot as plt


topo, lat, lon = rg.readasc(topo_file_name)
#topo, lat, lon = rg.readxyz('raw_topo/TOPO_MADEIRA_WGS84_UTM28_25m.dat')
topo[topo < 0] = 0
#plt.imshow(topo)

if (lon[0] >= 180) or (lat[0] >= 90):

    # subset
    north_lat, north_idx = nearest.nearest(lat, north)
    south_lat, south_idx = nearest.nearest(lat, south)
    east_lon, east_idx = nearest.nearest(lon, east)
    west_lon, weast_idx = nearest.nearest(lon, west)
    topo = topo[south_idx:north_idx, weast_idx:east_idx]
    lat, lon = lat[south_idx:north_idx], lon[weast_idx:east_idx]

elif (-180 <= lon[0] <= 180) or (-90 <= lat[0] <= 90):

    # subset
    north_lat, north_idx = nearest.nearest(lat, north)
    south_lat, south_idx = nearest.nearest(lat, south)
    east_lon, east_idx = nearest.nearest(lon, east)
    west_lon, weast_idx = nearest.nearest(lon, west)
    topo = topo[south_idx:north_idx, weast_idx:east_idx]
    lat, lon = lat[south_idx:north_idx], lon[weast_idx:east_idx]
#    lat_orig, lon_orig = lat, lon

    myProj = Proj("+proj=utm +zone=" + zone + ", +north +ellps=WGS84 +datum=WGS84 +units=m +no_defs")

    # inferior left and superior right points...
    lat_utm = []
    for lt in lat :
        lat_utm.append(myProj(lon[0], lt)[1])
#        lat_utm.append(transform(inProj, outProj, lon[0], lt)[1])
    lat_utm = np.array(lat_utm)

    lon_utm = []
    for lo in lon :
        lon_utm.append(myProj(lo, lat[0])[0])
#        lon_utm.append(transform(inProj, outProj, lo, lat[0])[0])
    lon_utm = np.array(lon_utm)

    lon = lon_utm
    lat = lat_utm


# change topo to resolution dx and dy
#topo_res = rgc.grid_points_change(topo, int(nx), int(ny), 'linear')
topo_res = rgc.grid_res_change(topo, lon, lat, dx, dy, 'nearest')
topo_res[topo_res < .5] = 0
lon = np.arange(lon[0], lon[-1], dx)
lat = np.arange(lat[0], lat[-1], dy)
#plt.imshow(topo_res)


x_dist = lon[-1]-lon[0] #cd.calc_distance_convert_utm(lat[0], lon[0], lat[0], lon[-1], zone = '28')
y_dist = lat[-1]-lat[0] #cd.calc_distance(lat[0], lon[0], lat[-1], lon[0])

nx = x_dist/dx
ny = y_dist/dy


#verifica se o numero é odd/impar e torna-o par para melhor divisao dos cores. Depois muda as dimensões
if not bool(int(nx-1) & 1) :
    nx = nx -1
    lon = lon[:-1]
    topo_res = topo_res[:, :-1]
    
if not bool(int(ny-1) & 1) :
    ny = ny -1
    lat = lat[:-1]
    topo_res = topo_res[:-1, :]
    
topo_res = dz * np.round(topo_res/dz)

plt.pcolormesh(topo_res);plt.colorbar(); plt.show()

print('Number o grid points x, y = ' + str(topo_res.shape[1]-1) + ', ' + str(topo_res.shape[0]-1))

nc_output = Dataset(case_name + '_static', 'w')
# Global Attributes
nc_output.description = 'Contains static information like orography'
nc_output.history = 'Created by Ricardo Faria - OOM (ricardo88faria@gmail.com) ' + time.ctime(time.time())
nc_output.source = 'netCDF4 python'
nc_output.source = 'netCDF4 python'
nc_output.origin_lat = np.float(-17) #lat[0]
nc_output.origin_lon = np.float(32.0) #lon[0]
nc_output.origin_z = np.float(0)
nc_output.origin_x = np.float(0)
nc_output.origin_y = np.float(0)
nc_output.rotation_angle = np.float(0)
nc_output.origin_time = '2011-06-16 00:00:00 +00'

nc_output.createDimension('x', nx)
nc_output.createDimension('y', ny)


nc_x = nc_output.createVariable('x',  np.float32, 'x')
nc_y = nc_output.createVariable('y',  np.float32, 'y')
nc_topo = nc_output.createVariable('zt',  np.float32, ('y', 'x'), fill_value = -9999.0, zlib=True)

nc_x.units = 'm'
nc_x.long_name = 'distance to origin in x-direction'
nc_x.axis = 'x'
nc_y.units = 'm'
nc_y.long_name = 'distance to origin in y-direction'
nc_y.axis = 'y'


nc_topo.long_name = 'terrain_height'
nc_topo.res_orig = np.float32(dx) #str(dx) + 'x' + str(dy) + ' m'
nc_topo.source = 'RAM'
nc_topo.units = 'm'
#nc_topo.coordinates = 'E_UTM N_UTM lon lat'
#nc_topo.grid_mapping = 'crsUTM: E_UTM N_UTM crsETRS: lon lat'

# set let bottom corner to 0 meters
lon = lon - lon[0]
lat = lat - lat[0]

nc_x[:] = np.arange(lon[1]/2, dx*nx+1, dx) # lon #np.arange(0, dx*nx, dx)
nc_y[:] = np.arange(lat[1]/2, dy*ny+1, dy) # lat #np.arange(0, dy*ny, dy)
nc_topo[:] = topo_res[:-1,:-1] 

nc_output.close()

print('Process finished! Topo is created with following specs in centered grids: \n [nx, ny, nz] = ' +  str([nx-1, ny-1, nz]) + ' \n [dx, dy, dz] = ' +  str([dx, dy, dz]))

#cfg = open(case_name + '.cfg', 'w')
#cfg.write(' nx = ' + str(nx-1) + '\n ny = ' + str(ny-1) + '\n nz = ' + str(nz))
#cfg.write('\n dx = ' + str(dx-1) + '\n dy = ' + str(dy-1) + '\n dz = ' + str(dz-1))
#cfg.close()

cfg = pd.DataFrame({'dx': [dx], 'dy': [dy], 'dz': [dz], 'nx': [nx], 'ny': [ny], 'nz': [nz], 'north': [north], 'south': [south], 'east': [east], 'west': [west]}) #, 'north_lat': [north_lat], 'south_lat': [south_lat], 'east_lon': [east_lon], 'west_lon': [west_lon]
cfg.to_csv(case_name + '.cfg', index=None)  # str(dx) + '_' + str(dy) +

