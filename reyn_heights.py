 #!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jan 29 08:18:07 2023

@author: sarahdennis
"""
import numpy as np
import random

from domain import Height

#------------------------------------------------------------------------------
# Example Height Functions
#------------------------------------------------------------------------------

class RandomHeight(Height):
    def __init__(self, x0, xf, N, h_min, h_max):
        h_str = "Discrete Height"

        
        Nx = (xf-x0)*N + 1
        hs = np.zeros(Nx)
        for i in range (Nx):
            # hs[i] = h_min + (h_max - h_min) * random.random()
            hs[i] = h_min + (h_max - h_min) * random.random()/(i+1)
        y0 = 0
        yf = max(hs)
            
        U = 1    
        h_avg = np.mean(hs)
        Re = U*h_avg
        super().__init__(x0, xf, y0, yf, Nx-1, hs, U, Re, h_str)
  

#------------------------------------------------------------------------------   

class SinsusoidalHeight(Height): 
    #h(x) = h_min + r(1 + cos(kx))
    def __init__(self, x0, xf, N, h_avg, r, k):
        Nx = (xf-x0)*N + 1
        self.h_mid = h_avg
        self.r = r 
        self.k = k
 
        h_str = "Sinusoidal Height"
            
        dx = (xf - x0)/(Nx-1)
        xs = np.asarray([x0 + i*dx for i in range(Nx)])
        hs = np.asarray([self.h_fun(x) for x in xs])
        
        y0 = 0
        yf = 1.1*(h_avg+r) 

        U = 1    
        Re = U*h_avg
        super().__init__(x0, xf, y0, yf, N, hs, U, Re, h_str)

    # def h_fun(self, x):
    #     return self.h_mid * (1 + self.r * np.cos(self.k*x))    
    
    # def h_fun(self, x):
    #     return self.h_mid * (1 + self.r * np.cos(self.k*x) + self.r/self.k * np.sin(self.k*2*x))    
    
    def h_fun(self, x):
        return self.h_mid * (1 + 4*self.r * np.cos(self.k*x) + self.r/2*self.k * np.sin(self.k*x))    

#------------------------------------------------------------------------------
class ConstantHeight(Height):
    def __init__(self, x0, xf, N, h0):
        h_str = "Constant Height"
        Nx = (xf-x0)*N + 1
        hs = np.ones(Nx)*h0
        
        y0 = 0
        yf = 1.1*h0
        U = 1    
        Re = U*h0
        super().__init__(x0, xf, y0, yf, N, hs, U, Re, h_str)

#------------------------------------------------------------------------------
class LinearHeight(Height): #slider bearing
    
    def __init__(self, x0, xf, N, h0, hf):
        self.h0 = h0
        self.hf = hf
        self.x0 = x0
        self.m = (hf - h0)/(xf - x0)
        self.N_regions=1
        h_str = "Wedge Slider"
        Nx = (xf-x0)*N + 1
        dx = (xf - x0)/(Nx-1)
        xs = np.asarray([x0 + i*dx for i in range(Nx)])
        hs = np.asarray([self.h_fun(x) for x in xs])

        y0 = 0
        yf = max(h0, hf)
        U = 1    
        h_avg = np.mean(hs)
        Re = U*h_avg
        super().__init__(x0, xf, y0, yf, N, hs, U, Re, h_str)

    def h_fun(self, x):
        return self.h0 + self.m * (x - self.x0)

#------------------------------------------------------------------------------    
class StepWaveHeight(Height): # uniform width [h1, h2, ..., hN+1]

    def __init__(self, x0, xf, N, N_steps, h_steps):
        self.N_steps = N_steps
        self.h_steps = h_steps
        self.step_width = (xf - x0)/(N_steps+1)
        h_str = "N=%d Step Height"%N_steps
        
        y0 = 0
        yf = max(h_steps)
        Nx = (xf-x0)*N + 1
        hs = self.make_hs(x0, xf, Nx, N_steps, h_steps, self.step_width)
        U = 1    
        h_avg = np.mean(hs)
        Re = U*h_avg
        super().__init__(x0, xf, y0, yf, N, hs, U, Re, h_str)

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
    def __init__(self, x0, xf, N, N_regions, x_peaks, h_peaks):
        self.h_peaks = h_peaks
        self.x_peaks = x_peaks
        self.N_regions = N_regions #=len(hpeaks)-1
        
        h_str = "Piecewise Linear"
        Nx = (xf-x0)*N + 1

        hs, self.slopes, self.widths = self.make_hs(x0, xf, Nx, x_peaks, h_peaks)
        y0 = 0
        yf = max(hs)
        U = 1    
        h_avg = np.mean(hs)
        Re = U*h_avg
        super().__init__(x0, xf, y0, yf, N, hs, U, Re, h_str)
        
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

    
#-----------------------------------------------------------------------------

class PiecewiseLinearHeight(Height):
    def __init__(self, x0, xf, N, N_regions, x_peaks, h_peaks):
        self.h_peaks = h_peaks
        self.x_peaks = x_peaks
        self.N_regions = N_regions #=len(hpeaks)-1
        
        h_str = "Piecewise Linear (new)"
        Nx = (xf-x0)*N + 1

        hs, self.slopes, self.widths = self.make_hs(x0, xf, Nx, x_peaks, h_peaks)
        y0 = 0
        yf = max(hs)
        U = 1     
        Re = U*np.mean(hs)
        super().__init__(x0, xf, y0, yf, N, hs, U, Re, h_str)
        
    def make_hs(self, x0, xf, Nx, x_peaks, h_peaks):
        slopes = np.zeros(self.N_regions)
        widths = np.zeros(self.N_regions)
        
        for r in range(self.N_regions):            

            slopes[r] = (h_peaks[r+1,0] - h_peaks[r,1])/(x_peaks[r+1] - x_peaks[r])

            
        hs = np.zeros(Nx)
        dx = (xf - x0)/(Nx-1)
        r = 0
        for i in range(Nx):
            xi = x0 + i*dx
            
            if xi > x_peaks[r+1] and r+1 < self.N_regions:
                r +=1
                
            widths[r] = xi - x_peaks[r]
            hs[i] = h_peaks[r,1] + slopes[r] * (xi - x_peaks[r])

        return  hs, slopes, widths

    
    
    
    
    
    
    
    
    
    