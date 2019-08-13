#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 22 11:31:46 2019

@author: ricardofaria

Calculate Geostrophic wind profile from wrf file area.
Based on:
- Practical_Meteorology-v1.02b-WholeBookColor pag.302 
- https://unidata.github.io/python-gallery/examples/Ageostrophic_Wind_Example.html

"""


__all__ = ['geostr', 'geostr_met']



def geostr_profile(theta_2d, height_cross, u_cross, v_cross) :
    
    '''

    Calculate Geostrophic wind map (2d) from 2d array area.
    https://unidata.github.io/python-gallery/examples/Ageostrophic_Wind_Example.html
    
    Parameters
    ----------
    array_2d_geoheight : read numpy array [m]
    array_1d_lat : read numpy array [deg]
    array_1d_lon : read numpy array [deg]

    Returns
    -------
    array_2d : return interplated and extrapolated 2d array

    '''
    
    
    import numpy as np
    import atmosphere as atm
    from nearest import nearest as near
    
    
    rib = atm.Rib(theta_2d[:], theta_2d[0], height_cross[:], height_cross[0], u_cross[:], v_cross[:])
    idx_stable_boundary_layer_height = near(rib[:], 0.25)[1]
    stable_boundary_layer_height = height_cross[idx_stable_boundary_layer_height]
    
    
    return(stable_boundary_layer_height, idx_stable_boundary_layer_height)



def geostr_met(array_2d_geoheight, array_1d_lat, array_1d_lon) :
    
    '''

    Calculate Geostrophic wind map (2d) from 2d array area.
    https://unidata.github.io/python-gallery/examples/Ageostrophic_Wind_Example.html
    
    Parameters
    ----------
    array_2d_geoheight : read numpy array [m]
    array_1d_lat : read numpy array [deg]
    array_1d_lon : read numpy array [deg]

    Returns
    -------
    array_2d : return interplated and extrapolated 2d array

    '''
    
    
    import numpy as np
    import metpy.calc as mpcalc
    from metpy.units import units
    
    
    # Combine 1D latitude and longitudes into a 2D grid of locations
    lon_2d, lat_2d = np.meshgrid(array_1d_lon, array_1d_lat)
    
    # Set up some constants based on our projection, including the Coriolis parameter and
    # grid spacing, converting lon/lat spacing to Cartesian
    f = mpcalc.coriolis_parameter(np.deg2rad(lat_2d)).to('1/s')
    dx, dy = mpcalc.lat_lon_grid_deltas(lon_2d + 360, lat_2d)
    dx, dy = np.array(dx), np.array(dy)
    dy *= -1
    
    # In MetPy 0.5, geostrophic_wind() assumes the order of the dimensions is (X, Y),
    # so we need to transpose from the input data, which are ordered lat (y), lon (x).
    # Once we get the components,transpose again so they match our original data.
    geo_wind_u, geo_wind_v = mpcalc.geostrophic_wind(array_2d_geoheight.data * units.m, f, dx, dy)
    
    return(geo_wind_u, geo_wind_v)
    

def geostr(array_2d_press, array_2d_temp, array_1d_lat, array_1d_lon) :
    
    '''

    Calculate Geostrophic wind profile (1 point value representing input 2d array area).
    Based on Practical_Meteorology-v1.02b-WholeBookColor pag.302
    
    Parameters
    ----------
    array_2d_press : read numpy array [Pa]
    array_2d_temp : read numpy array [k]
    array_1d_lat : read numpy array [deg]
    array_1d_lon : read numpy array [deg]

    Returns
    -------
    array : return interplated and extrapolated value

    '''
    
    
    import numpy as np
    import metpy.calc as mpcalc
    import geopy.distance
    
    
    # Combine 1D latitude and longitudes into a 2D grid of locations
    lon_2d, lat_2d = np.meshgrid(array_1d_lon[[0, -1]], array_1d_lat[[0, -1]])
    x_grd_center = np.int(array_2d_press.shape[1]/2)
    y_grd_center = np.int(array_2d_press.shape[0]/2)
    
    # Set up some constants based on our projection, including the Coriolis parameter and
    # grid spacing, converting lon/lat spacing to Cartesian
    f = (1.4584e-4)*np.sin(np.deg2rad(array_1d_lat[y_grd_center]))
#    f = mpcalc.coriolis_parameter(np.deg2rad(lat_2d))
#    dx, dy = mpcalc.lat_lon_grid_deltas(lon_2d, lat_2d)
#    dy *= -1
    dx = geopy.distance.vincenty((array_1d_lat[0], array_1d_lon[x_grd_center]), (array_1d_lat[-1], array_1d_lon[x_grd_center])).m
    dy = geopy.distance.vincenty((array_1d_lat[y_grd_center], array_1d_lon[0]), (array_1d_lat[y_grd_center], array_1d_lon[-1])).m
    
    rho_mean = np.nanmean(rho(array_2d_temp, array_2d_press))
    
    geo_wind_u = (-1 / (rho_mean * f)) * ((np.nanmean(array_2d_press[-1, :]) - np.nanmean(array_2d_press[0, :])) / dy)
    geo_wind_v = (1 / (rho_mean * f)) * ((np.nanmean(array_2d_press[:, -1]) - np.nanmean(array_2d_press[:, 0])) / dx)
    
    geo_wind = np.array([geo_wind_u, geo_wind_v])
    
#    # Combine 1D latitude and longitudes into a 2D grid of locations
#    lon_2d, lat_2d = np.meshgrid(array_1d_lon, array_1d_lat)
#    x_grd_center = np.int(array_2d_press.shape[1]/2)
#    y_grd_center = np.int(array_2d_press.shape[0]/2)
#    
#    f = (1.4584e-4)*np.sin(np.deg2rad(array_1d_lat[y_grd_center]))   # np.sin(array_1d_lat[y_grd_center]*np.pi/180)
#    
#    dx = geopy.distance.vincenty((array_1d_lat[0], array_1d_lon[x_grd_center]), (array_1d_lat[1], array_1d_lon[x_grd_center])).m
#    dy = geopy.distance.vincenty((array_1d_lat[y_grd_center], array_1d_lon[0]), (array_1d_lat[y_grd_center], array_1d_lon[1])).m
#
#    rho_mean = rho(array_2d_temp, array_2d_press)
#    
#    geo_wind_u = np.zeros((rho_mean.shape[0]-1, rho_mean.shape[1]-1))
#    geo_wind_v = np.zeros((rho_mean.shape[0]-1, rho_mean.shape[1]-1))
#    for x in range(rho_mean.shape[1]-1) :
#        for y in range(rho_mean.shape[0]-1) :
##            print(x, y)
#            geo_wind_u[y, x] = (-1 / (rho_mean[y, x] * f)) * ((array_2d_press[y+1, x] - array_2d_press[y, x]) / dy)
#            geo_wind_v[y, x] = (1 / (rho_mean[y, x] * f)) * ((array_2d_press[y, x+1] - array_2d_press[y, x]) / dx)
    
#    np.nanmean(geo_wind_u), np.nanmean(geo_wind_v)
#    plt.pcolormesh(geo_wind_u); plt.colorbar(); plt.show()
#    plt.pcolormesh(geo_wind_v); plt.colorbar(); plt.show()
    
    return(geo_wind)
    

def rho(T, p):
    
    """
    
    Calculates air density (rho)
    
    """
    
    
    Rd = 287.

#    Tv   = T * (1+0.61*qv) # Virtual temperature

    rho = p / (Rd * T) # Air density [kg m^-3]

    return rho

