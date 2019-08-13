#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Created on Sat Aug 25 11:27:33 2018

@author: ricardofaria

Script to create static driver netcdf input for PALM following PALM Input Data Standard (PIDS) v1.9
Edit my_test_setup_static.nc – contains all static information like orography, buildings, and set surface classification.

"""


landuse_file_name = 'raw_landuse/clc_06.tif'
case_name = 'madeira_breeze_50m_offline_nesting'
zone = '28'
interp_mode = 'linear'  # linear
nsurface_fraction = 3
#zlad = 8

building_height = 10
#canopy_height = 5
#alpha_lad = 3
#beta_lad = 3


#import matplotlib.pyplot as plt
import copy
import numpy as np
from netCDF4 import Dataset #, num2date
#import glob
import pandas as pd
from pyproj import Proj, transform
import read_geo as rg
#import calc_distance as cd
import res_grid_change as rgc
import nearest
import CLClv5_to_PALM


cfg = pd.read_csv(case_name + '.cfg')
dx = cfg.dx.values[0]
dy = cfg.dy.values[0]
dz = cfg.dz.values[0]
nx = cfg.nx.values[0]
ny = cfg.ny.values[0]
zlad = nz = cfg.nz.values[0]
north = cfg.north.values[0]
south = cfg.south.values[0]
east = cfg.east.values[0]
west = cfg.west.values[0]
north_lat = cfg.north_lat.values[0]
south_lat = cfg.south_lat.values[0]
east_lon = cfg.east_lon.values[0]
west_lon = cfg.west_lon.values[0]

vegetation_type_specs = pd.read_csv('vegetation_type_specs.csv', sep=';')


lu, lat, lon = rg.readgeotiff(landuse_file_name)
lu = lu[::-1, :]
lu[lu <= 0] = 52311


if (-180 <= lon[0] <= 180) or (-90 <= lat[0] <= 90):

    # subset
    north_idx = nearest.nearest(lat, north)[1]
    south_idx = nearest.nearest(lat, south)[1]
    east_idx = nearest.nearest(lon, east)[1]
    weast_idx = nearest.nearest(lon, west)[1]
    lu = lu[south_idx:north_idx, weast_idx:east_idx]
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


#lu_res = rgc.grid_res_change(lu, lon, lat, dx, dy, 'nearest')
#plt.pcolormesh(lu_res), plt.colorbar()

x_new_points = nx # np.floor((lon[-1]-lon[0])/dx) #lon.shape[0]/(dx/(lon[1]-lon[0]))
y_new_points = ny # np.floor((lat[-1]-lat[0])/dy) #lat.shape[0]/(dx/(lat[1]-lat[0]))
lu_res = rgc.grid_points_change(lu, x_new_points+1, y_new_points+1, 'nearest')
#plt.pcolormesh(lu_res), plt.colorbar()

if lu_res.shape[1]-1 != nx :
    lu_res = lu_res[:, 0:lu_res.shape[1]-1]
if lu_res.shape[0]-1 != ny :
    lu_res = lu_res[0:lu_res.shape[0]-1, :]
print('Number o grid points x, y = ' + str(lu_res.shape[1]-1) + ', ' + str(lu_res.shape[0]-1))

lon = np.arange(lon[0], lon[-1], dx)
lat = np.arange(lat[0], lat[-1], dy)

x_dist = lon[-1]-lon[0] #cd.calc_distance_convert_utm(lat[0], lon[0], lat[0], lon[-1], zone = '28')
y_dist = lat[-1]-lat[0] #cd.calc_distance(lat[0], lon[0], lat[-1], lon[0])

nx = x_dist/dx
ny = y_dist/dy


vegetation = np.zeros_like(lu_res)
vegetation[:] = lu_res[:]
vegetation = CLClv5_to_PALM.reclassify(vegetation, 'vegetation')
#plt.pcolormesh(vegetation, vmin =0); plt.colorbar()

pavement = np.zeros_like(lu_res)
pavement[:] = lu_res[:]
pavement = CLClv5_to_PALM.reclassify(pavement, 'pavement')

building = np.zeros_like(lu_res)
building[:] = lu_res[:]
building = CLClv5_to_PALM.reclassify(building, 'building')

water = np.zeros_like(lu_res)
water[:] = lu_res[:]
water = CLClv5_to_PALM.reclassify(water, 'water')
#water[water == 0] = 3

soil = np.zeros_like(lu_res)
soil[:] = -9999.0
soil[vegetation > 0] = 2
soil[pavement > 0] = 2
#soil[:] = lu_res[:]
#soil = CLC_to_PALM.reclassify(soil, 'soil')

#surface_fraction = np.zeros_like(lu_res)
surface_fraction = np.array([np.zeros_like(lu_res), np.zeros_like(lu_res), np.zeros_like(lu_res)])
surface_fraction[:] = -9999.0
for i in range(surface_fraction.shape[0]):
    surface_fraction[i][vegetation > 0] = 0
    surface_fraction[i][pavement > 0] = 0
    surface_fraction[i][water > 0] = 0

surface_fraction[0][vegetation > 0] = 1
surface_fraction[1][pavement > 0] = 1
surface_fraction[2][water > 0] = 1
#surface_fraction[:] = lu_res[:]
#surface_fraction = CLC_to_PALM.reclassify(surface_fraction, 'surface_fraction')

buildings_2d = np.zeros_like(lu_res)
buildings_2d[:] = -9999.0
buildings_2d[building > 0] = building_height # 30 metros de altura ao solo


building_id = np.zeros_like(lu_res)
building_id[:] = building
#id_nr = np.arange(1, len(np.where(building_id > 0)[0]))
count = 0
for x in range(building_id.shape[1]) :
    for y in range(building_id.shape[0]) :
        if building_id[y, x] > 0 :
            
            count = count + 1
            building_id[y, x] = count

# LAD profile following Markkanen et al. (2003)
z_range_lad = np.append(0, np.arange(dz/2, np.nanmax(vegetation_type_specs.iloc[:, 13])+1, dz))
lad = np.zeros(shape=(int(len(z_range_lad)), lu_res.shape[0], lu_res.shape[1]))
lad[:] = -9999.0
#lad[0][vegetation == 4] = 0.1

z_range_buildings_3d = np.append(0, np.arange(dz/2, building_height+1, dz))
buildings_3d = np.zeros(shape=(int(len(z_range_buildings_3d)), lu_res.shape[0], lu_res.shape[1]))

# calc 3d buildings and lad profile for each x and y points (3d)
for x in range(lad.shape[2]) :
    for y in range(lad.shape[1]) :
        
        # lad profiles
        if vegetation[y, x] > 0 :
            
            veg_id = vegetation[y, x]
            alpha_lad = vegetation_type_specs.iloc[int(veg_id-1), 14]
            beta_lad = vegetation_type_specs.iloc[int(veg_id-1), 15]
            canopy_height = vegetation_type_specs.iloc[int(veg_id-1), 13]
            lai = vegetation_type_specs.iloc[int(veg_id-1), 2]
            
            z_range_lad_xy = np.append(0, np.arange(dz/2, canopy_height+1, dz))
            
            for k_idx, k in enumerate(z_range_lad_xy) :
                
                int_bpdf = (( z_range_lad[k_idx] / canopy_height )**(alpha_lad - 1))*(( 1 - ( z_range_lad[k_idx] / canopy_height ) )**(beta_lad-1)) * (( z_range_lad[k_idx]+dz-z_range_lad[k_idx] ) / canopy_height)
                pre_lad = lai * ( ( ( z_range_lad[k_idx] / canopy_height )**( alpha_lad-1 ) ) * ( ( 1 - ( z_range_lad[k_idx] / canopy_height ) )**( beta_lad-1 ) ) / int_bpdf ) / canopy_height
#                lai_beta * ( ( ( zw(k) / canopy_height )**( alpha_lad-1.0_wp ) ) * ( ( 1.0_wp - ( zw(k) / canopy_height ) )**( beta_lad-1.0_wp ) )
                lad[int(k_idx), y, x] = pre_lad #0.5 * ( pre_lad(k-1) + pre_lad(k) )
            
        
        # 3d buildings fill buildings height with 1
        if buildings_2d[y, x] > 0 :
            
#            building_height = buildings_2d[y, x]
            z_range_buildings_3d_xy = np.append(0, np.arange(dz/2, building_height+1, dz))
            
            for k_idx, k in enumerate(z_range_buildings_3d_xy) :
                
                buildings_3d[int(k_idx), y, x] = 1
        


nc_output = Dataset(case_name + '_static', 'r+')

try:
    nc_output.createDimension('nsurface_fraction', nsurface_fraction)
    nc_output.createDimension('zlad', len(z_range_lad))
    nc_output.createDimension('z', len(z_range_buildings_3d))
    
    nc_nsurface_fraction = nc_output.createVariable('nsurface_fraction', np.int32, 'nsurface_fraction')
    nc_zlad = nc_output.createVariable('zlad', np.float32, 'zlad')
    nc_z = nc_output.createVariable('z', np.float32, 'z')
    
    nc_lad = nc_output.createVariable('lad',  np.float32, ('zlad', 'y', 'x'), fill_value = -9999.0, zlib=True)
    nc_vegetation = nc_output.createVariable('vegetation_type',  np.byte, ('y', 'x'), fill_value = -127.0, zlib=True)
    nc_pavement = nc_output.createVariable('pavement_type',  np.byte, ('y', 'x'), fill_value = -127.0, zlib=True)
    nc_building = nc_output.createVariable('building_type',  np.int32, ('y', 'x'), fill_value = -127.0, zlib=True)
    nc_water = nc_output.createVariable('water_type',  np.byte, ('y', 'x'), fill_value = -127.0, zlib=True)
    nc_soil = nc_output.createVariable('soil_type',  np.byte, ('y', 'x'), fill_value = -127.0, zlib=True)
    nc_surface_fraction = nc_output.createVariable('surface_fraction',  np.float32, ('nsurface_fraction', 'y', 'x'), fill_value = -9999.0, zlib=True)
    nc_buildings_2d = nc_output.createVariable('buildings_2d',  np.float32, ('y', 'x'), fill_value = -9999.0, zlib=True)
    nc_building_id = nc_output.createVariable('building_id',  np.float32, ('y', 'x'), fill_value = -9999, zlib=True)
    nc_buildings_3d = nc_output.createVariable('buildings_3d',  np.float32, ('z', 'y', 'x'), fill_value = -127.0, zlib=True)
    
    nc_z.units = 'm'
    nc_z.long_name = 'height above origin'
    #nc_z.standard_name = 'height_above_mean_sea_level'
    nc_z.positive = 'up'
    nc_z.axis = 'z'

    nc_nsurface_fraction.units = ''
    nc_zlad.units = 'm'
    
    nc_lad.long_name = 'leaf_area_density'
    nc_lad.res_orig = np.float32(dx) #str(dx) + 'x' + str(dy) + ' m'
    nc_lad.source = 'RAM'
    nc_lad.units = 'm2 m-3'
    
    nc_vegetation.long_name = 'vegetation_type_classification'
    #nc_vegetation.res_orig = np.float32(dx) #str(dx) + 'x' + str(dy) + ' m'
    nc_vegetation.source = 'RAM'
    nc_vegetation.units = ''
    
    nc_pavement.long_name = 'pavement_type_classification'
    #nc_pavement.res_orig = np.float32(dx) #str(dx) + 'x' + str(dy) + ' m'
    nc_pavement.source = 'RAM'
    nc_pavement.units = ''
    
    nc_building.long_name = 'building_type_classification'
    nc_building.res_orig = np.float32(dx) #str(dx) + 'x' + str(dy) + ' m'
    nc_building.source = 'RAM'
    nc_building.units = ''
    
    nc_water.long_name = 'water_type_classification'
    #nc_water.res_orig = np.float32(dx) #str(dx) + 'x' + str(dy) + ' m'
    nc_water.source = 'RAM'
    nc_water.units = ''
    
    nc_soil.long_name = 'soil_type_classification'
    #nc_soil.res_orig = np.float32(dx) #str(dx) + 'x' + str(dy) + ' m'
    nc_soil.source = 'RAM'
    nc_soil.units = ''
    #nc_soil.lod = np.float32(1)
    
    nc_surface_fraction.long_name = 'surface_fraction'
    #nc_surface_fraction.res_orig = np.float32(dx) #str(dx) + 'x' + str(dy) + ' m'
    nc_surface_fraction.source = 'RAM'
    nc_surface_fraction.units = ''
    
    nc_buildings_2d.long_name = 'building_height'
    nc_buildings_2d.res_orig = np.float32(dx) #str(dx) + 'x' + str(dy) + ' m'
    nc_buildings_2d.lod = np.float32(1)
    nc_buildings_2d.source = 'RAM'
    nc_buildings_2d.units = 'm'
    
    nc_building_id.long_name = 'building_id_numbers'
    nc_building_id.res_orig = np.float32(dx) #str(dx) + 'x' + str(dy) + ' m'
    nc_building_id.source = 'RAM'
    nc_building_id.units = ''
    
    nc_buildings_3d.long_name = 'building_height'
    nc_buildings_3d.res_orig = np.float32(dx)
    nc_buildings_3d.lod = np.float32(2)
    nc_buildings_3d.source = 'RAM'
    nc_buildings_3d.units = '1'
    

except:
    nc_nsurface_fraction = nc_output.variables['nsurface_fraction']
    nc_zlad = nc_output.variables['zlad']
    nc_z = nc_output.variables['z']
    nc_lad = nc_output.variables['lad']
    nc_vegetation = nc_output.variables['vegetation_type']
    nc_pavement = nc_output.variables['pavement_type']
    nc_building = nc_output.variables['building_type']
    nc_water = nc_output.variables['water_type']
    nc_soil = nc_output.variables['soil_type']
    nc_surface_fraction = nc_output.variables['surface_fraction']
    nc_buildings_2d = nc_output.variables['buildings_2d']
    nc_building_id = nc_output.variables['building_id']
    nc_buildings_3d = nc_output.variables['buildings_3d']


vegetation[vegetation == -9999.0] = -127.0
vegetation[vegetation == 0] = -127.0
pavement[pavement == -9999.0] = -127.0
pavement[pavement == 0] = -127.0
building[building == -9999.0] = -127.0
building[building == 0] = -127.0
water[water == -9999.0] = -127.0
water[water == 0] = -127.0
soil[soil == -9999.0] = -127.0
soil[soil == 0] = -127.0
buildings_2d[buildings_2d == -9999.0] = -9999.0
buildings_2d[buildings_2d == 0] = -9999.0
building_id[building_id == -9999.0] = -9999
building_id[building_id == 0] = -9999

#surface_fraction[surface_fraction == -9999.0] = -9999
#surface_fraction[surface_fraction == 0] = -9999
#plt.pcolormesh(building)


vegetation_tmp = copy.copy(vegetation[:])
vegetation_tmp[vegetation_tmp ==-127.0] = 0
pavement_tmp = copy.copy(pavement[:])
pavement_tmp[pavement_tmp ==-127.0] = 0
building_tmp = copy.copy(building[:])
building_tmp[building_tmp ==-127.0] = 0
water_tmp = copy.copy(water[:])
water_tmp[water_tmp ==-127.0] = 0
soil_tmp = copy.copy(soil[:])
soil_tmp[soil_tmp ==-127.0] = 0
total = vegetation_tmp+pavement_tmp+building_tmp+water_tmp+soil_tmp

land = vegetation_tmp + pavement_tmp + building_tmp

nc_nsurface_fraction[:] = np.arange(0, nsurface_fraction)
nc_zlad[:] = z_range_lad #np.arange(dz/2, dz*nz, dz)
nc_z[:] = z_range_buildings_3d

nc_lad[:] = lad[:,:-1,:-1]
nc_vegetation[:] = vegetation[:-1,:-1]
nc_pavement[:] = pavement[:-1,:-1]
nc_building[:] = building[:-1,:-1]
nc_water[:] = water[:-1,:-1]
nc_soil[:] = soil[:-1,:-1]
nc_surface_fraction[:] = surface_fraction[:,:-1,:-1]
nc_buildings_2d[:] = buildings_2d[:-1,:-1]
nc_building_id[:] = building_id[:-1,:-1]
nc_buildings_3d[:] = buildings_3d[:,:-1,:-1]

topo = nc_output.variables['zt'][:]
topo[water[:-1,:-1] == 3] = 0
nc_output.variables['zt'][:] = topo

nc_output.close()

print('Process finished! LandUse is introduced with following specs in centered grids: \n [nx, ny, nz] = ' +  str([nx-1, ny-1, nz]) + ' \n [dx, dy, dz] = ' +  str([dx, dy, dz]))
