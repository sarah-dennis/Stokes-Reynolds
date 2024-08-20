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

    def __init__(self, x0, xf, y0, yf, N):
        self.x0 = x0            # left boundary x = x0
        self.xf = xf            # right boundary x = xf
        self.y0 = y0            # lower surface y = y0 
        self.yf = yf            # upper surface y = height(x) <= yf

        self.N = N  # scale: number of grid points in [0,1]
        self.Nx = int((xf-x0)*N + 1) #total number of x grid points
        self.Ny = int((yf-y0)*N + 1) #total number of y grid points

        self.dx = 1/N
        self.dy = 1/N
        
        self.xs = np.linspace(x0, xf, self.Nx)
        self.ys = np.linspace(y0, yf, self.Ny)
        
#------------------------------------------------------------------------------
# Domain for Reynolds solver
class Height(Domain):

    def __init__(self, x0, xf, y0, yf, N, hs, U, Re, h_str):
        super().__init__(x0, xf, y0, yf, N) # -> {dx, dy, xs, ys}
        
        self.hs = hs
        self.hxs = center_diff(self.hs, self.Nx, self.dx)
        self.hxxs = center_second_diff(self.hs, self.Nx, self.dx)
        
        self.h_str = h_str
        self.h_eq = "h(x)"

        self.h_max = max(self.hs)
        self.h_min = min(self.hs)
        
        self.U = U   # flat boundary velocity

        self.visc = 1 # dynamic viscosity 

        self.Re = Re # Re = U h / L
 
# Domain for Stokes solver
class Space(Domain):
    def __init__(self, x0, xf, y0, yf, N, U, Q, Re, h_str):
        super().__init__(x0, xf, y0, yf, N)
        self.filestr=h_str
        self.U = U   # flat boundary velocity
        self.flux = Q
        self.visc = 1 # dynamic viscosity 

        self.Re = Re # Re = U h / L
        
    def set_space(self, grid):
        self.space = grid
        
#------------------------------------------------------------------------------
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
   
    
    
    