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
# PWL Height
#------------------------------------------------------------------------------
class PWL_Height(Height):
    def __init__(self, x0, xf, N, N_regions, x_peaks, h_peaks, U, dP, filestr):
        self.h_peaks = h_peaks
        self.x_peaks = x_peaks
        self.N_regions = N_regions # = len(h_peaks)-1
        
        hs, self.slopes, self.widths = self.make_hs(x0, xf, N, x_peaks, h_peaks)
        
        y0 = 0
        yf = max(hs)  
        
        super().__init__(x0, xf, y0, yf, N, hs, U, dP, filestr)
        
    def make_hs(self, x0, xf, N, x_peaks, h_peaks):
        slopes = np.zeros(self.N_regions)
        widths = np.zeros(self.N_regions)
        
        for r in range(self.N_regions):            

            slopes[r] = (h_peaks[r+1,0] - h_peaks[r,1])/(x_peaks[r+1] - x_peaks[r])

        Nx = (xf-x0)*N + 1
        hs = np.zeros(Nx)
        dx = 1/N
        r = 0
        for i in range(Nx):
            xi = x0 + i*dx
            
            if xi > x_peaks[r+1] and r+1 < self.N_regions:
                r +=1
                
            widths[r] = xi - x_peaks[r]
            hs[i] = h_peaks[r,1] + slopes[r] * (xi - x_peaks[r])
        return  hs, slopes, widths
    

    
    
#------------------------------------------------------------------------------
# Other Height Functions
#------------------------------------------------------------------------------

class RandomHeight(Height):
    def __init__(self, x0, xf, N, h_min, h_max, U, dP, filestr):
        # h_str = "./examples/" +f"Rand_H{h_max}_U{U}_dP{dP}_N{N}"

        Nx = (xf-x0)*N + 1
        hs = np.zeros(Nx)
        for i in range (Nx):
            # hs[i] = h_min + (h_max - h_min) * random.random()
            hs[i] = h_min + (h_max - h_min) * random.random()/(i+1)
        y0 = 0
        yf = max(hs)
        
        super().__init__(x0, xf, y0, yf, Nx-1, hs, U, dP, filestr)
  

#------------------------------------------------------------------------------   
class SinusoidalHeight(Height): 
    #h(x) = h_min + r(1 + cos(kx))
    def __init__(self, x0, xf, N, h_avg, r, k, U, dP, filestr):
        Nx = (xf-x0)*N + 1
        self.h_mid = h_avg
        self.r = r 
        self.k = k
 
        # h_str = "./examples/" + f"sin_h{h_avg}_r{r}_k{k}_U{U}_dP{dP}_N{N}"
            
        dx = 1/N
        xs = np.asarray([x0 + i*dx for i in range(Nx)])
        hs = np.asarray([self.h_fun(x) for x in xs])
        
        y0 = 0
        yf = (h_avg+r) 

        super().__init__(x0, xf, y0, yf, N, hs, U, dP, filestr)

    def h_fun(self, x):
        return self.h_mid * (1 + self.r * np.cos(self.k*x))   
    
#------------------------------------------------------------------------------    
class CircleHeight(Height):
    
    def __init__(self, x0, xf, N, r, h0, l, U, dP, filestr):
        Nx = (xf-x0)*N + 1
        dx = 1/N
        
        self.h0 = h0
        self.r = r 
        self.l = l
    
        
            
        dx = (xf - x0)/(Nx-1)
        xs = np.asarray([x0 + i*dx for i in range(Nx)])
        hs = np.asarray([self.h_fun(x) for x in xs])
        
        y0 = 0
        yf = h0 + r
        # h_str = "./examples/" +f"circ_r{r}_H{yf}_l{l}_U{U}_dP{dP}_N{N}"
        super().__init__(x0, xf, y0, yf, N, hs, U, dP, filestr)

    def h_fun(self, x):
        # return self.h0 + np.sqrt(self.r**2 - (x-self.r)**2)
        if x <= self.l:
            return self.h0+self.r
        elif x <= self.l + 2*self.r:
            return self.h0 + self.r - np.sqrt(self.r**2 - (x-self.l-self.r)**2)
        else:
            return self.h0+self.r
    
    
#------------------------------------------------------------------------------
class ConstantHeight(Height):
    def __init__(self, x0, xf, N, h0, filestr):
        # h_str = "./examples/" +f"Cnsnt_H{h0}_L{xf-x0}_N{N}"
        Nx = (xf-x0)*N + 1
        hs = np.ones(Nx)*h0
        
        y0 = 0
        yf = 1.1*h0
        U = 1    

        dP=2
        super().__init__(x0, xf, y0, yf, N, hs, U, dP, filestr)

#------------------------------------------------------------------------------
class StepHeight(Height):
    def __init__(self, x0, xf, N, h0, hf, x_step, U,dP, filestr):
        self.x_step = x_step
        self.N_steps = 1
        self.h_steps= [h0, hf]
        self.step_width = (xf - x0)/2
        # h_str = "./examples/" +f"BFS_H{hf}_L{xf-x0}_U{U}_dP{dP}_N{N}"
        y0 = 0
        yf = max(self.h_steps)
        Nx = (xf-x0)*N + 1
        hs = self.make_hs(x0, xf, Nx, self.N_steps, self.h_steps, self.step_width)  


        super().__init__(x0, xf, y0, yf, N, hs, U, dP, filestr)

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

    
    
    
    
    
    
    