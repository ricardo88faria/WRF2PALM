#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Created on Thu Nov 30 12:16:19 2017

@author: ricardofaria

Read geo *.asc file and retreive data 2d with lat lon 1d

"""


import numpy as np
import linecache
import pandas as pd
#import math
#import matplotlib.pyplot as plt

__all__ = ["readasc", "readxyz", "readgeotiff"]


def get_num(x):

    return float(''.join(ele for ele in x if ele.isdigit() or ele == '.' or ele == '-'))


def readasc(file) :

    """

    read a geo .asc file from file handle

    Parameters
    ----------
    file : file handle
           file to read from

    Returns
    -------
    xyz : namedtuple
        returns a named tuple with data and lat lon coords.

    """

#    file = 'srtm_33_06.asc'

    data_array = np.loadtxt(file, skiprows=6)
    data_array = data_array[::-1, :]

    ncols = linecache.getline(file, 1)
    ncols = get_num(ncols)
    nrows = linecache.getline(file, 2)
    nrows = get_num(nrows)
    xllcorner = linecache.getline(file, 3)
    xllcorner = get_num(xllcorner)
    yllcorner = linecache.getline(file, 4)
    yllcorner = get_num(yllcorner)
    cellsize = linecache.getline(file, 5)
    cellsize = get_num(cellsize)
    NODATA_value = linecache.getline(file, 6)
    NODATA_value = get_num(NODATA_value)

    lat = np.arange(yllcorner, yllcorner+(cellsize*nrows), cellsize)
    lon = np.arange(xllcorner, xllcorner+(cellsize*ncols), cellsize)

#    plt.pcolormesh(data_array)

    return(data_array, lat, lon)


def readxyz(file) :

    """

    read a geo .dat file from file handle

    Parameters
    ----------
    file : file handle
           file to read from

    Returns
    -------
    xyz : namedtuple
        returns a named tuple with data and lat lon coords.
    To implement, creation of a empty array with xyz mat size and fill each y and x point, the rest can be set to nan
    """

    #file = 'lat_lon_corvo50m.dat'

    data_array = pd.read_csv(file, sep = ' ', header = None)
    lon, lat = np.unique(data_array[0]), np.unique(data_array[1])

#    data_array = np.meshgrid(data_array[0], data_array[1], data_array[2])
#    data_array = np.vstack(np.meshgrid(data_array[0], data_array[1], data_array[2])).reshape(3,-1).T
#    yx = zip(lat, lon)
#    data_array_list = data_array.values.tolist()

    data_array['idx'] = data_array.index

#    n2 = len(data_array_list)
#    n = int(math.sqrt(n2))
#    ndarray = np.zeros([n, n])
#    ndarray_1d_idx = np.arange(len(lat), len(lat)*len(lat)+1, len(lat))

#    ndarray = np.arange(len(data_array_list))
    
#    ndarray = np.reshape(data_array[2][::-1], (len(lat), len(lon)))
    ndarray = np.array(data_array[2][::-1]).reshape((len(lat), len(lon)))
    
    data_array = ndarray[::-1, ::-1]

#    plt.pcolormesh(data_array)

#    count_rows = 0
#    count_cols = 0
#    for idx, l in enumerate(data_array_list) :  # enumerate(yx)
##        print(idx)
##        ndarray[l] = data[idx]
#        if idx in ndarray_1d_idx :
#
#            ndarray[count_rows, count_cols-1] = data_array[2][idx]
#            count_rows = count_rows+1
#            count_cols = 0
#            print(count_cols, count_rows, idx)
#
#        elif idx not in ndarray_1d_idx :
#
#            ndarray[count_rows, count_cols-1] = data_array[2][idx]
#            count_cols = count_cols+1
#            print(count_cols, count_rows, idx)

    return(data_array, lat, lon)


def readgeotiff(file) :

    """

    read a geo .tif file from file handle

    Parameters
    ----------
    file : file handle
           file to read from

    Returns
    -------
    xyz : namedtuple
        returns a named tuple with data and lat lon coords.

    """

#    import rasterio
#
##    file = 'raw_topo/srtm_33_06.tif'
#    r = rasterio.open(str(file))
#
#    lat = np.linspace(r.bounds.bottom, r.bounds.top, r.height)
#    lon = np.linspace(r.bounds.left, r.bounds.right, r.width)
#
#    r = r.affine
#    data_array = r.read()[0,]
#
##    plt.pcolormesh(data_array)
#
#    return(data_array, lat, lon)
    
    from osgeo import gdal
    
    ds = gdal.Open(str(file), gdal.GA_ReadOnly)
#    rb = ds.GetRasterBand(0)
    
    width = ds.RasterXSize
    height = ds.RasterYSize
    
#    source_layer = ds.GetLayer()
    pos = ds.GetGeoTransform()
    x_min = pos[0]
    y_min = pos[3] + width*pos[4] + height*pos[5] 
    x_max = pos[0] + width*pos[1] + height*pos[2]
    y_max = pos[3] 
    
    lat = np.linspace(y_min, y_max, height)
    lon = np.linspace(x_min, x_max, width)
    
    data_array = ds.ReadAsArray()
    
    return(data_array, lat, lon)
    

    