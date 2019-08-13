#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 10 10:52:23 2018

@author: ricardofaria

Interpolate and extrapolate 2d data with nan's using:
    * Radial Basis Function Interpolation / Kernel Smoothing;
    * Gaussian Process Regression / Kriging.

"""


__all__ = ['kriging', 'kern_smooth', 'interp3d']



def interp3d(vals, var, arr) :

    import numpy as np

    '''

    Interpolate 3d array to vals levels in arr.
    var (ex: wind vel) in function to arr (ex: height) and mask values outside max height

    '''


    n_var = np.empty((len(vals), var.shape[-2], var.shape[-1]))
#    depth_ratio = np.linspace(0, 1, depth.shape[0])

    vals = np.array(vals)
    for x in range(var.shape[-1]) :
        for y in range(var.shape[-2]) :

            var_2d = np.array(var[:, y, x])
            arr_2d = np.array(arr[:, y, x])

            n_var[:, y, x] = np.interp(vals, arr_2d, var_2d)
            mask_vals_idx = np.where(vals > arr_2d.max())
            n_var[mask_vals_idx, y, x] = np.nan

    return(n_var)


def kriging(array_in, x, y, theta0 = 0.1, thetaL = .001, thetaU = 1., nugget = 0.01) :

    """

    read numpy array with nans and return interplated and extrapolated 2d array using Gaussian Process Regression / Kriging method

    Parameters
    ----------
    array_in : read numpy array with nans

    Returns
    -------
    array_out : return interplated and extrapolated 2d array

    """

    import numpy as np
#    from sklearn.gaussian_process import GaussianProcess
    from sklearn.gaussian_process import GaussianProcessRegressor
    from sklearn.gaussian_process.kernels import RBF, ConstantKernel as C
    
    kernel = C(1.0, (1e-3, 1e3)) * RBF([5,5], (1e-2, 1e2))
    gp = GaussianProcessRegressor(kernel=kernel, n_restarts_optimizer=1)
#    gp = GaussianProcess(theta0=theta0, thetaL=thetaL, thetaU=thetaU, nugget=nugget)

    xx, yy = np.meshgrid(x, y)
    vals = ~np.isnan(array_in)

    gp.fit(X=np.column_stack([xx[vals],yy[vals]]), y=array_in[vals])

    xx_yy_as_cols = np.column_stack([xx.flatten(), yy.flatten()])

    array_out = gp.predict(xx_yy_as_cols).reshape(array_in.shape)

#    plt.imshow(GD1,interpolation='nearest')

    return(array_out)


def kern_smooth(array_in, x, y, method = 'linear') :

    """

    read numpy array with nans and return interplated and extrapolated 2d array using Radial Basis Function Interpolation / Kernel Smoothing.

    Parameters
    ----------
    array_in : read numpy array with nans

    Returns
    -------
    array_out : return interplated and extrapolated 2d array

    """

    import numpy as np
    import scipy.interpolate as interpolate

    xx, yy = np.meshgrid(x, y)
    vals = ~np.isnan(array_in)

    f = interpolate.Rbf(xx[vals], yy[vals], array_in[vals])
    array_out = f(xx, yy)

#    plt.imshow(GD1,interpolation='nearest')

    return(array_out)


##tests
#xx, yy = np.meshgrid(range(left_u[t, :].shape[1]), range(left_u[t, :].shape[0]))
#vals = ~np.isnan(left_u[t, :])
#array_in = left_u[t, :]
#
#%%timeit
#f = interpolate.griddata(xx[vals], yy[vals], array_in[vals], method='linear')
#array_out = f(xx, yy)
#
#points = np.vstack((xx.ravel(), yy.ravel())).T
#values = array_in.ravel()
#%%timeit
#rbfi = interpolate.Rbf(points[:,0],points[:,1],values)
#dl = rbfi(xx.ravel(),yy.ravel())
