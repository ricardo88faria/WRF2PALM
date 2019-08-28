#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 26 10:33:40 2019

@author: ricardofaria

numpy operations

"""



def using_clump(a):
    
    
    """
    
    split:
        [nan,nan, 1 , 2 , 3 , nan, nan, 10, 11 , nan, nan, nan, 23, 1, nan, 7, 8]
    
    Into:
        [[1,2,3], [10,11], [23,1], [7,8]]
    
    """
    
    import numpy as np
    
#    slices = np.ma.clump_unmasked(np.ma.masked_invalid(a))
#    
#    for s in slices :
#        b = a[s]
        
    vals = [a[s] for s in np.ma.clump_unmasked(np.ma.masked_invalid(a))]
#    idx = [s for s in np.ma.clump_unmasked(np.ma.masked_invalid(a))]
    
    return vals

