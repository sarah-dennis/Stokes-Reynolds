#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 25 17:33:35 2023

@author: sarahdennis
"""
import numpy as np
#------------------------------------------------------------------------------

# Domain: 2D grid on [x0, xf] and [y0, yf] on Nx and Ny grid points

class Domain:

    def __init__(self, x0, xf, y0, yf, Nx, Ny):
        self.x0 = x0            # left boundary x = x0
        self.xf = xf            # right boundary x = xf
        self.y0 = y0             # lower surface y = y0 
        self.yf = yf            # upper surface y = height(x) <= yf
        
        self.Nx = Nx  #Number of x grid points
        self.Ny = Ny  #Number of y grid points
        
        self.dx = (self.xf - self.x0)/(self.Nx-1)
        self.dy = (self.yf - self.y0)/(self.Ny-1)
        
        self.xs = np.asarray([self.x0 + i*self.dx for i in range(self.Nx)])
        self.ys = np.asarray([self.y0 + i*self.dy for i in range(self.Ny)])
             
    # def get_index_x(self, x):
    #     return int((x-self.x0)//self.dx)
    
    
def center_diff(fs, N, dx):
    D_lower = -1*np.ones(N)
    D_upper = np.ones(N)
    D = np.diagflat(D_lower[1:N], -1) + np.diagflat(D_upper[0:N-1], 1)
    
    D = D/(2*dx)
    fs_dx = D@fs 
    
    fs_dx[0] = fs_dx[1]
    fs_dx[-1] = fs_dx[-2]
    
    return np.asarray(fs_dx)
      
def center_second_diff(fs, N, dx):
    D_lower = np.ones(N-1)
    D_upper = np.ones(N-1)
    D_center = -2*np.ones(N)
    D = np.diagflat(D_lower, -1) +  np.diagflat(D_center, 0)+ np.diagflat(D_upper, 1)
    
    D[0][N-1] = 0
    D[N-1][0] = 0
        
    D = D/(dx**2)
    fs_dxx = D@fs 
    
    return np.asarray(fs_dxx)
   
    
    
    