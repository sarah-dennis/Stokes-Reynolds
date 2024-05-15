 #!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jan 29 08:18:07 2023

@author: sarahdennis
"""
import numpy as np
import random

import domain as dm
#------------------------------------------------------------------------------
# A domain [x0, xf],[y0, yf] with a height function y = h_fun(x) 
class Height(dm.Domain):

    def __init__(self, x0, xf, y0, yf, Nx, Ny, hs, h_str, h_eq):

        super().__init__(x0, xf, y0, yf, Nx, Ny) # -> {dx, dy, xs, ys}
        
        self.hs = hs
        self.hxs = dm.center_diff(self.hs, self.Nx, self.dx)
        self.hxs = dm.center_second_diff(self.hs, self.Nx, self.dx)
        
        self.h_str = h_str
        self.h_eq = h_eq

        self.h_max = max(self.hs)
        self.h_min = min(self.hs)
        self.h_avg = np.mean(self.hs)
        
        self.U = 1   # velocity on to flat lower surface
        self.dens = 1 # density
        self.visc = 1 # dynamic viscosity 
    
        self.Re = self.U * self.h_avg / self.visc
 

#------------------------------------------------------------------------------
# Example Height Functions
#------------------------------------------------------------------------------

class RandomHeight(Height):
    def __init__(self, x0, xf, Nx, h_min, h_max):
        h_str = "Discrete Height"
        h_eq = "h(x) = h_i"
        
        hs = np.zeros(Nx)
        for i in range (Nx):
            hs[i] = h_min + (h_max - h_min) * random.random()
            y0 = 0
            yf = 1.1*max(hs)
            Ny = Nx
        super().__init__(x0, xf, y0, yf, Nx, Ny, hs, h_str, h_eq)
  

#------------------------------------------------------------------------------   

class SinsusoidalHeight(Height): 
    #h(x) = h_min + r(1 + cos(kx))
    def __init__(self, x0, xf, Nx, h_avg, r, k):
        
        self.h_mid = h_avg
        self.r = r 
        self.k = k

        h_eq = "h(x) = %0.1f + %0.1f(1 + \cos(%d x))"%(h_avg, r, k) 
        h_str = "Sinusoidal Height"
            
        dx = (xf - x0)/(Nx-1)
        xs = np.asarray([x0 + i*dx for i in range(Nx)])
        hs = np.asarray([self.h_fun(x) for x in xs])
        
        y0 = 0
        yf = 1.1*(h_avg+r) 
        Ny = Nx
        
        super().__init__(x0, xf, y0, yf, Nx, Ny, hs, h_str, h_eq)

    def h_fun(self, x):
        return self.h_mid * (1 + self.r * np.cos(self.k*x))    
#------------------------------------------------------------------------------
class ConstantHeight(Height):
    def __init__(self, x0, xf, Nx, h0):
        h_str = "Constant Height"
        h_eq = "h(x) = %.2f"%h0
        
        hs = np.ones(Nx)*h0
        
        y0 = 0
        yf = 1.1*h0
        Ny = Nx
        
        super().__init__(x0, xf, y0, yf, Nx, Ny, hs, h_str, h_eq)

#------------------------------------------------------------------------------
class LinearHeight(Height): #slider bearing
    
    def __init__(self, x0, xf, Nx, h0, hf):
        self.h0 = h0
        self.hf = hf
        self.x0 = x0
        self.m = (hf - h0)/(xf - x0)
        
        h_eq = "h(x) = %0.1f + %0.1f(x - %0.1f)"%(self.h0, self.m, self.x0)
        h_str = "Wedge Slider"
        
        dx = (xf - x0)/(Nx-1)
        xs = np.asarray([x0 + i*dx for i in range(Nx)])
        hs = np.asarray([self.h_fun(x) for x in xs])

        y0 = 0
        yf = 1.1*max(h0, hf)
        Ny = Nx
        
        super().__init__(x0, xf, y0, yf, Nx, Ny, hs,  h_str, h_eq)

    def h_fun(self, x):
        return self.h0 + self.m * (x - self.x0)

#------------------------------------------------------------------------------

#------------------------------------------------------------------------------    
class StepWaveHeight(Height): # uniform width [h1, h2, ..., hN+1]

    def __init__(self, x0, xf, Nx, N_steps, h_steps):
        self.N_steps = N_steps
        self.h_steps = h_steps
        self.step_width = (xf - x0)/(N_steps+1)
        h_str = "N=%d Step Height"%N_steps
        h_eq = "h(x)"
        
        hs = self.make_hs(x0, xf, Nx, N_steps, h_steps, self.step_width)
        
        y0 = 0
        yf = 1.1*max(h_steps)
        Ny = Nx 

        super().__init__(x0, xf, y0, yf, Nx, Ny, hs, h_str, h_eq)

    def make_hs(self, x0, xf, Nx, n_steps, h_steps, step_width):
        hs = np.zeros(Nx)
        index_width = Nx / (n_steps + 1)

        j=0
        for i in range(Nx):

            if i >= (j+1)*index_width :
                
                j += 1
            hs[i] = h_steps[j]
 
        return hs

class StepHeight(StepWaveHeight):
    def __init__(self, x0, xf, Nx, h0, hf, x_step):
        self.x_step = x_step
        N_steps = 1
        h_steps = [h0, hf]
        super().__init__(x0, xf, Nx, N_steps, h_steps)
#-----------------------------------------------------------------------------
    
class SawtoothHeight(Height):
    def __init__(self, x0, xf, Nx, N_regions, x_peaks, h_peaks):
        self.h_peaks = h_peaks
        self.x_peaks = x_peaks
        self.N_regions = N_regions #=len(hpeaks)-1
        hs, self.slopes, self.widths = self.make_hs(x0, xf, Nx, x_peaks, h_peaks)
                
        h_eq = "h(x)"
        h_str = "Piecewise Linear"
        
        y0 = 0
        yf = 1.1*max(hs)
        Ny = Nx 
        
        super().__init__(x0, xf, y0, yf, Nx, Ny, hs,  h_str, h_eq)
        
    def make_hs(self, x0, xf, Nx, x_peaks, h_peaks):
        slopes = np.zeros(self.N_regions)
        widths = np.zeros(self.N_regions)
        
        for r in range(self.N_regions):            

            slopes[r] = (h_peaks[r+1] - h_peaks[r])/(x_peaks[r+1] - x_peaks[r])

            
        hs = np.zeros(Nx)
        dx = (xf - x0)/(Nx-1)
        r = 0
        for i in range(Nx):
            xi = x0 + i*dx
            
            if xi > x_peaks[r+1] and r+1 < self.N_regions:
                r +=1
                
            widths[r] = xi - x_peaks[r]
            hs[i] = h_peaks[r] + slopes[r] * (xi - x_peaks[r])

        return  hs, slopes, widths
 

    
    
    
    
    
    
    
    
    
    