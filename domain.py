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
        
        if BC == "periodic": #periodic
            self.dx = (xf - x0)/(Nx)
        elif BC == "fixed": #fixed
            self.dx = (xf - x0)/(Nx-1)
        
        self.xs = [x0 + i*self.dx for i in range(Nx)]
  
    
def center_diff(fs, domain):
    D_lower = -1*np.ones(domain.Nx)
    D_upper = np.ones(domain.Nx)
    D = np.diagflat(D_lower[1:domain.Nx], -1) + np.diagflat(D_upper[0:domain.Nx-1], 1)
    
    if domain.BC == "periodic": #periodic 
        D[0][domain.Nx-1] = -1
        D[domain.Nx-1][0] = 1
        
    elif domain.BC == "fixed": #prescribed 
        None # boundary heights are not used
        
    D = D/(2*domain.dx)
        
    fs_dx = D@fs 
 
    return fs_dx
      
def center_second_diff(fs, domain):
    D_lower = np.ones(domain.Nx-1)
    D_upper = np.ones(domain.Nx-1)
    D_center = -2*np.ones(domain.Nx)
    D = np.diagflat(D_lower, -1) +  np.diagflat(D_center, 0)+ np.diagflat(D_upper, 1)
    
    if domain.BC == 0: #periodic 
        D[0][domain.Nx-1] = 1
        D[domain.Nx-1][0] = 1
        
    elif domain.BC == 1: #prescribed 
        #None # boundary heights are not used
        D[0][domain.Nx-1] = 0
        D[domain.Nx-1][0] = 0
        
    D = D/(domain.dx**2)
        
    fs_dxx = D@fs 
    
    return fs_dxx
   
    

    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    