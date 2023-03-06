#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 25 17:33:35 2023

@author: sarahdennis
"""
import numpy as np

class Domain:

    def __init__(self, x0, xf, eta, U, Nx, BC):
        self.x0 = x0
        self.xf = xf
        self.Nx = Nx
        self.BC = BC
        self.eta = eta
        self.U = U
        
        if BC == 0: #periodic
            self.dx = (xf - x0)/(Nx)
        elif BC == 1: #fixed
            self.dx = (xf - x0)/(Nx-1)
        
        self.xs = [x0 + i*self.dx for i in range(Nx)]
  
    
def center_diff(fs, Domain):
    D_lower = -1*np.ones(Domain.Nx)
    D_upper = np.ones(Domain.Nx)
    D = np.diagflat(D_lower[1:Domain.Nx], -1) + np.diagflat(D_upper[0:Domain.Nx-1], 1)
    
    if Domain.BC == 0: #periodic 
        D[0][Domain.Nx-1] = -1
        D[Domain.Nx-1][0] = 1
        
    elif Domain.BC == 1: #prescribed 
        None # boundary heights are not used
        
    D = D/(2*Domain.dx)
        
    fs_dx = D@fs 
 
    return fs_dx
      
def center_second_diff(fs, Domain):
    D_lower = np.ones(Domain.Nx-1)
    D_upper = np.ones(Domain.Nx-1)
    D_center = -2*np.ones(Domain.Nx)
    D = np.diagflat(D_lower, -1) +  np.diagflat(D_center, 0)+ np.diagflat(D_upper, 1)
    
    if Domain.BC == 0: #periodic 
        D[0][Domain.Nx-1] = 1
        D[Domain.Nx-1][0] = 1
        
    elif Domain.BC == 1: #prescribed 
        #None # boundary heights are not used
        D[0][Domain.Nx-1] = 0
        D[Domain.Nx-1][0] = 0
        
    D = D/(Domain.dx**2)
        
    fs_dxx = D@fs 
 
    #graph.plot_2D_multi([fs_dxx[1:-1], fs[1:-1]], xs[1:-1], "discretized", ["f_xx", "f"])
    return fs_dxx
   
    

    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    