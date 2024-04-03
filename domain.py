#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 25 17:33:35 2023

@author: sarahdennis
"""
import numpy as np

class Domain:
##TODO switch/convert eta to Re
    def __init__(self, x0, xf, U, Nx, Ny, BC='fixed'):
        self.x0 = x0
        self.xf = xf
        self.Nx = Nx
        self.Ny = Ny
        self.BC = BC
        self.eta = 1 #dynamic viscosity in reynolds eq
        self.U = U
        
        if BC == "periodic": #periodic
            self.dx = (xf - x0)/(Nx)
        else: #fixed
            self.dx = (xf - x0)/(Nx-1)
        
        self.xs = np.asarray([x0 + i*self.dx for i in range(Nx)])


 # #------------------------------------------------------------------------------
    # # Reynolds number
    # Re = U * height.h_avg / visc
    # # print("Re = %.4f"%Re)
    # #TODO: separate domain instances for X and Y
    
    def set_ys(self, height, Ny):
        # always fixed BC in y
        # used by Velocity once height is made
        self.Ny = Ny
        self.yf = height.h_max
        self.y0 = height.h_min
        self.dy = (self.yf - self.y0)/(Ny-1)
        self.ys = np.asarray([self.y0 + i*self.dy for i in range(Ny)])
  
    def get_index(self, x):
        return int((x-self.x0)//self.dx)
    
def center_diff(fs, domain):
    D_lower = -1*np.ones(domain.Nx)
    D_upper = np.ones(domain.Nx)
    D = np.diagflat(D_lower[1:domain.Nx], -1) + np.diagflat(D_upper[0:domain.Nx-1], 1)
    
    if domain.BC == "periodic": #periodic 
        D[0][domain.Nx-1] = -1
        D[domain.Nx-1][0] = 1
        
    D = D/(2*domain.dx)
    fs_dx = D@fs 
    
    if  domain.BC == "fixed": #prescribed 
        fs_dx[0] = fs_dx[1]
        fs_dx[-1] = fs_dx[-2]
    
    return np.asarray(fs_dx)
      
def center_second_diff(fs, domain):
    D_lower = np.ones(domain.Nx-1)
    D_upper = np.ones(domain.Nx-1)
    D_center = -2*np.ones(domain.Nx)
    D = np.diagflat(D_lower, -1) +  np.diagflat(D_center, 0)+ np.diagflat(D_upper, 1)
    
    if domain.BC == "periodic": #periodic 
        D[0][domain.Nx-1] = 1
        D[domain.Nx-1][0] = 1
        
    elif domain.BC == "fixed": #prescribed 
        # None # boundary not used
        D[0][domain.Nx-1] = 0
        D[domain.Nx-1][0] = 0
        
    D = D/(domain.dx**2)
    fs_dxx = D@fs 
    
    return np.asarray(fs_dxx)
   
    
    
    