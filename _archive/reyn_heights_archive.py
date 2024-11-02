# -*- coding: utf-8 -*-
"""
Created on Tue Aug  6 09:06:16 2024

@author: sarah
"""
import numpy as np

#------------------------------------------------------------------------------    
class StepWaveHeight(Height): # uniform width [h1, h2, ..., hN+1]

    def __init__(self, x0, xf, N, N_steps, h_steps, U):
        self.N_steps = N_steps
        self.h_steps = h_steps
        self.step_width = (xf - x0)/(N_steps+1)
        h_str = "N=%d Step Height"%N_steps
        
        y0 = 0
        yf = max(h_steps)
        Nx = (xf-x0)*N + 1
        hs = self.make_hs(x0, xf, Nx, N_steps, h_steps, self.step_width)  
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

    