# -*- coding: utf-8 -*-
"""
Created on Wed Aug 30 12:20:23 2023

@author: sarah
"""
import numpy as np
import domain as dm
import boundary as bc

# 2D Reynolds velicity field from 1D pressure (dp/dx -> Q )
def make_reyn_velocity(height, BC, pressure):

    u = np.zeros((height.Ny, height.Nx))
    v = np.zeros((height.Ny, height.Nx))
    
    U = BC.U
    if isinstance(BC, bc.Fixed):
       Q = get_reyn_flux(height, BC, pressure)
        
    elif isinstance(BC, bc.Mixed):
        Q = BC.Q 
        
    for i in range(height.Nx):

        h = height.hs[i]
        hx = height.hxs[i]

        for j in range(height.Ny):
            y = height.ys[j]
            if y <= height.hs[i]:
                u[j,i] = (h-y)*(U*(h-3*y)/h**2 + 6*Q*y/h**3)
                v[j,i] = -2*hx * y**2 * (h-y) *(U/h**3 - 3*Q/h**4)
                
            else:
                u[j,i] = 0
                v[j,i] = 0
                
    return u, v

def get_reyn_flux(height, BC, pressure):
    h = height.hs[1]
    px = dm.center_first(height.dx, pressure.ps_1D[0:3])
    U = BC.U
    Q = (U*h)/2 - (px *(h**3))/12 #/visc
    return Q
