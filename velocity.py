# -*- coding: utf-8 -*-
"""
Created on Wed Aug 30 12:20:23 2023

@author: sarah
"""
import numpy as np

def getFluidVelocity(domain, height, pressure):

    # direction components
    v_x = np.zeros(domain.Nx)
    v_y = np.zeros(domain.Nx) 

    h = height.h_max
    for i in range(domain.Nx):
        y = height.ys[i]
        
        v_x[i] = 1/(2*domain.eta) * pressure.pxs[i] * (y**2 - y*h) + y*domain.U/h
        v_y[i] = domain.U * height.hxs[i]

    return v_x, v_y