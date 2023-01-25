#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 25 17:45:17 2023

@author: sarahdennis
"""
import numpy as np

  
def center_diff(f, Domain):
    D_lower = -1*np.ones(Domain.Nx)
    D_upper = np.ones(Domain.Nx)
    D = np.diagflat(D_lower[1:Domain.Nx], -1) + np.diagflat(D_upper[0:Domain.Nx-1], 1)
    
    if Domain.BC == 0: #periodic 
        D[0][Domain.Nx-1] = -1
        D[Domain.Nx-1][0] = 1
        
    elif Domain.BC == 1: #prescribed 
        None # boundary heights are not used
        
    D = D/(2*Domain.dx)
    fs = [f(x) for x in Domain.xs]
        
    fs_dx = D@fs 
 
    #graph.plot_2D_multi([fs_dx[1:-1], fs[1:-1]], xs[1:-1], "discretized", ["f_x", "f"])
    return fs_dx 
      
def center_second_diff(f, Domain):
    D_lower = np.ones(Domain.Nx-1)
    D_upper = np.ones(Domain.Nx-1)
    D_center = -2*np.ones(Domain.Nx)
    D = np.diagflat(D_lower, -1) +  np.diagflat(D_center, 0)+ np.diagflat(D_upper, 1)
    
    if Domain.BC == 0: #periodic 
        D[0][Domain.Nx-1] = 1
        D[Domain.Nx-1][0] = 1
        
    elif Domain.BC == 1: #prescribed 
        None # boundary heights are not used
        
    D = D/(Domain.dx**2)
    fs = [f(x) for x in Domain.xs]
        
    fs_dxx = D@fs 
 
    #graph.plot_2D_multi([fs_dxx[1:-1], fs[1:-1]], xs[1:-1], "discretized", ["f_xx", "f"])
    return fs_dxx
   
