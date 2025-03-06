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

    def __init__(self, x0, xf, y0, yf, N, dirstr=None, namestr=None):
        self.x0 = x0            # left boundary x = x0
        self.xf = xf            # right boundary x = xf 
        self.y0 = y0            # lower surface y = y0
        self.yf = yf            # upper boundary y = max h = yf-y0

        self.N = N  # scale: number of grid points in [0,1]
        self.Nx = int((xf-x0)*N + 1) #total number of x grid points
        self.Ny = int((yf-y0)*N + 1) #total number of y grid points
        
        self.dx = 1/N
        self.dy = 1/N
        
        self.xs = np.linspace(x0, xf, self.Nx)
        self.ys = np.linspace(y0, yf, self.Ny)
        self.dirstr = dirstr
        self.filestr= f"{dirstr}/{namestr}_N{N}"
        

#------------------------------------------------------------------------------
# Domain for Reynolds solver
class Height(Domain):

    def __init__(self, x0, xf, y0, yf, N, hs, U, dP, namestr):
        dirstr = f"./reyn_examples/{namestr}"
        super().__init__(x0, xf, y0, yf, N, dirstr,namestr) # -> {dx, dy, xs, ys}
        
        self.hs = hs
        self.hxs = center_diff(self.hs, self.Nx, self.dx)

        self.h_max = max(self.hs)
        self.h_min = min(self.hs)
        
        self.U = U    # velocity at flat boundary
        self.visc = 1#2.45/8.48 # kinematic viscosity (m^2/s)
        # Tavakol: density 848 kgm^-3, dynamic viscosity 2.45 Pas, flux 14.4mm^3/s
        self.p_ambient = 0 
        self.p0 = -dP
        self.pN = self.p_ambient
     
        self.Re = 0


    
# Domain for Stokes solver
class Space(Domain):
    def __init__(self, x0, xf, y0, yf, N, U, flux, Re, namestr):
        dirstr = f"./examples/{namestr}"
        super().__init__(x0, xf, y0, yf, N, dirstr, namestr)
        self.U = U    # velocity at flat boundary 
        self.visc = 1  # dynamic viscosity (Ns/m^2)
        self.dens = 1  # density (kg/m^3)
        # kinematic viscosity nu = dynamic viscosity mu / density rho
        
        self.p_ambient = 0 #10^5 Pa   
        self.flux=flux
        self.Re = Re #

#------------------------------------------------------------------------------
def center_diff(fs, Nx, dx):
    D_lower = -1*np.ones(Nx)
    D_upper = np.ones(Nx)
    D = np.diagflat(D_lower[1:Nx], -1) + np.diagflat(D_upper[0:Nx-1], 1)
    D = D/(2*dx)
    
    fs_dx = D@fs 
    
    # fs_dx[0] = fs_dx[1]
    fs_dx[0]= (-3*fs[0]+4*fs[1]-fs[2])/(2*dx)
    
    # fs_dx[-1] = fs_dx[-2]
    fs_dx[Nx-1]= (3*fs[Nx-1]-4*fs[Nx-2]+fs[Nx-3])/(2*dx)
    return np.asarray(fs_dx)
      
def center_second_diff(fs, Nx, dx):
    D_lower = np.ones(Nx-1)
    D_upper = np.ones(Nx-1)
    D_center = -2*np.ones(Nx)
    D = np.diagflat(D_lower, -1) +  np.diagflat(D_center, 0)+ np.diagflat(D_upper, 1)

    D = D/(dx**2)
    fs_dxx = D@fs 
    
    fs_dxx[0]=(2*fs[0]-5*fs[1]+4*fs[2]-fs[3])/(dx**2)
    fs_dxx[Nx-1]=(2*fs[Nx-1]-5*fs[Nx-2]+4*fs[Nx-3]-fs[Nx-4])/(dx**2)
    
    return np.asarray(fs_dxx)
   
    
    
    
