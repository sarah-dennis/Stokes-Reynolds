# -*- coding: utf-8 -*-
"""
Created on Wed Aug 28 11:58:25 2024

@author: sarah
"""
import numpy as np
import domain as dm

def vorticity(ex, u_2D, v_2D):
    uy_2D = np.gradient(u_2D, ex.dy, axis=0)
    vx_2D = np.gradient(v_2D, ex.dx, axis=1)
    w_2D = np.zeros((ex.Ny,ex.Nx))
    for j in range(ex.Ny):
        for i in range(ex.Nx):   
            w_2D[j,i] = vx_2D[j,i] - uy_2D[j,i]
    return w_2D

def pressure_gradient(ex, u_2D, v_2D):
    #TODO add BCs based on ex.space
    # eta coefficient 
    uy = np.gradient(u_2D, ex.dy, axis=0)
    uyy = np.gradient(uy, ex.dy, axis=0)
    ux = np.gradient(u_2D, ex.dx, axis=1)
    uxx = np.gradient(ux, ex.dx, axis=1)

    vy = np.gradient(v_2D, ex.dy, axis=0)
    vyy = np.gradient(vy, ex.dy, axis=0)
    vx = np.gradient(v_2D, ex.dx, axis=1)
    vxx = np.gradient(vx, ex.dx, axis=1)
    
    px =(uxx + uyy)
    py =(vxx + vyy)
    # print(px[:,0])
    # print(px[:,-1])
    # print(ex.dp_in)
    # print(ex.dp_out)
    return px, py

# def pressure(ex, px, py):
    