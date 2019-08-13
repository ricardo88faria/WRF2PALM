#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 20 12:32:22 2018

@author: ricardofaria

Reclassification of a 2d array Corine Land Cover land use categories level 5 to PALM following PALM Input Data Standard (PIDS) v1.9
Note that this setting is developed for mesoscale parametrization.

"""


def reclassify(array, classification) :


    '''

    reclassify(array, classification)
    return array with replaced values for palm, depending on classification that can be:
        'vegetation', 'pavement', 'building', 'water', 'soil' and 'surface_fraction'

    '''
    
    import pandas as pd
    import numpy as np
    
    tab = pd.read_csv('raw_landuse/COS_2_PALM_num.csv')
    
    if classification == 'vegetation' :
        for l in range(tab.shape[0]) :
            array[array == tab.iloc[l,0]] = tab.iloc[l,1]
        array[np.isnan(array)] = -9999.0
    elif classification == 'pavement' :
        for l in range(tab.shape[0]) :
            array[array == tab.iloc[l,0]] = tab.iloc[l,2]
        array[np.isnan(array)] = -9999.0
    elif classification == 'building' :
        for l in range(tab.shape[0]) :
            array[array == tab.iloc[l,0]] = tab.iloc[l,3]
        array[np.isnan(array)] = -9999.0
    elif classification == 'water' :
        for l in range(tab.shape[0]) :
            array[array == tab.iloc[l,0]] = tab.iloc[l,4]
        array[np.isnan(array)] = -9999.0
    elif classification == 'soil' :
        for l in range(tab.shape[0]) :
            array[array == tab.iloc[l,0]] = tab.iloc[l,5]
        array[np.isnan(array)] = -9999.0

    return(array)
