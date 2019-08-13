#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul  4 14:35:52 2019

@author: ricardofaria

Check PALM divisible number of processors for the number of grids

"""


nr_grids = 618*1218           # model grid number
mx_nr_pr = 1000             # max number of processors


print(nr_grids, 'number of grids are divisible by: \n')

for n in range(1, mx_nr_pr+1) :
    if nr_grids % n == 0 :
        print(n)
