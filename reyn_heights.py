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
    def __init__(self, x0, xf, N, N_regions, x_peaks, h_peaks, filestr):
        self.h_peaks = h_peaks
        self.x_peaks = x_peaks
        self.N_regions = N_regions # = len(h_peaks)-1
        
        hs, self.slopes, self.widths, i_peaks = self.make_hs(x0, xf, N, x_peaks, h_peaks)
        
        y0 = 0
        yf = max(hs)  

        super().__init__(x0, xf, y0, yf, N, hs, i_peaks, filestr)
        
    def make_hs(self, x0, xf, N, x_peaks, h_peaks):
        slopes = np.zeros(self.N_regions)
        widths = np.zeros(self.N_regions)
        i_peaks = np.zeros(self.N_regions+1, dtype=int)
        for r in range(self.N_regions):            

            slopes[r] = (h_peaks[r+1,0] - h_peaks[r,1])/(x_peaks[r+1] - x_peaks[r])

        Nx = int((xf-x0)*N + 1)
        hs = np.zeros(Nx)
        dx = 1/N
        r = 0
        for i in range(Nx):
            xi = x0 + i*dx
            
            if xi > x_peaks[r+1] and r+1 < self.N_regions:
                r +=1
            
            i_peaks[r+1] = i    
            widths[r] = xi - x_peaks[r]
            hs[i] = h_peaks[r,1] + slopes[r] * (xi - x_peaks[r])

        return  hs, slopes, widths, i_peaks
    
#------------------------------------------------------------------------------
# PWC Height
#------------------------------------------------------------------------------
class PWC_Height(PWL_Height):
    def __init__(self, x0, xf, N, N_regions, x_peaks, h_peaks, filestr):
        
        #solver requires minimum 3 regions
        while N_regions <3 :
            x_peak_mid = x_peaks[-1] - x_peaks[-2]/2
            x_peak_end = x_peaks[-1]
            x_peaks[-1] = x_peak_mid
            x_peaks = np.append(x_peaks, x_peak_end)
            h_peaks = np.append(h_peaks, h_peaks[-1])
            h_peaks = np.reshape(h_peaks, (N_regions+2,2))
            N_regions +=1
            
            
        self.h_steps = np.zeros(N_regions)
        for i in range(N_regions):
            self.h_steps[i] = h_peaks[i,1]

        super().__init__(x0, xf, N, N_regions, x_peaks, h_peaks, filestr)
        


#------------------------------------------------------------------------------
# Other Height Functions
#------------------------------------------------------------------------------

class RandomHeight(Height):
    def __init__(self, x0, xf, N, h_min, h_max, filestr):

        Nx = (xf-x0)*N + 1
        hs = np.zeros(Nx)
        for i in range (Nx):
            hs[i] = h_min + (h_max - h_min) * random.random()/(i+1)
        y0 = 0
        yf = max(hs)
        i_peaks = np.asarray(range(Nx))
        super().__init__(x0, xf, y0, yf, Nx-1, hs, i_peaks,filestr)
  

#------------------------------------------------------------------------------   
class SinusoidalHeight(Height): 
    #h(x) = h_min + r(1 + cos(kx))
    def __init__(self, x0, xf, N, H, h, filestr):
        Nx = (xf-x0)*N + 1
        self.H = H
        self.h = h
        self.L = (xf-x0)
        
        # h_str = "./examples/" + f"sin_h{h_avg}_r{r}_k{k}_U{U}_dP{dP}_N{N}"
        y0 = 0
        yf =max(h, H+h)
        
        dx = 1/N
        xs = np.asarray([x0 + i*dx for i in range(Nx)])
        hs = np.asarray([self.h_fun(x) for x in xs])
        
        i_peaks = [0, Nx-1]

        super().__init__(x0, xf, y0, yf, N, hs, i_peaks, filestr)

    def h_fun(self, x):
        return self.H/2 * (1 + np.cos(2*np.pi/self.L*(self.L/2 - x))) + self.h
    
class BumpHeight(Height): 
    #h(x) = h_min + r(1 + cos(kx))
    def __init__(self, x0, xf, N, lam, H, h0, filestr):
        Nx = (xf-x0)*N + 1
        # h_str = "./examples/" + f"sin_h{h_avg}_r{r}_k{k}_U{U}_dP{dP}_N{N}"
        y0 = 0
        
        self.H=H         
        self.x_scale = (xf-x0)/2  
        self.h0 = h0
        self.lam=lam
        
        dx = 1/N
        xs = np.asarray([x0 + i*dx for i in range(Nx)])
        hs = np.asarray([self.h_fun(x) for x in xs])
        yf = np.max(hs)
        i_peaks = [0, Nx-1]
       
        super().__init__(x0, xf, y0, yf, N, hs, i_peaks, filestr)
  
    def h_fun(self, x):
        return self.H*(1-(self.lam/2)*(1+np.cos(np.pi*x/self.x_scale)))-(self.H-self.h0)
   
class LogisticHeight(Height):
    def __init__(self, x0, xf, N, H, h, center, delta, filestr):
       
        Nx = (xf-x0)*N + 1
        dx = 1/N
        self.h = h
        self.H = H
        self.delta = delta #slope = delta*(H-h)/4
        self.center = center
        xs = np.asarray([x0 + i*dx for i in range(Nx)])
        hs = np.asarray([self.h_fun(x) for x in xs])  
        y0 = 0
        yf = np.max(hs)
        i_peaks = [0, Nx-1]
       
        super().__init__(x0, xf, y0, yf, N, hs, i_peaks, filestr)

    def h_fun(self, x):
        
        return (self.h + (self.H-self.h) / ( 1 + np.exp(self.delta*(x-self.center))))  

#------------------------------------------------------------------------------    
class CircleHeight(Height):
    
    def __init__(self, x0, xf, N, r, dxdr, h0, filestr):
        Nx = int((xf-x0)*N + 1)
        dx = 1/N
        
        self.h0 = h0 #clearance
        self.r = r 
        self.l = ((xf-x0)-2*r)/2
            
        self.dxdr =  dxdr
        
            
        dx = (xf - x0)/(Nx-1)

        xs = np.asarray([x0 + i*dx for i in range(Nx)])
        hs = np.asarray([self.h_fun(x) for x in xs])
        y0 = 0
        yf = max(hs)
        # print(self.l, self.l*N)
        if self.l != 0:
            i_peaks = np.asarray([0, (self.l+self.dxdr)*N, Nx-1-(self.l+self.dxdr)*N, Nx-1], int)
        else:
            i_peaks = [0, Nx-1]
       
        super().__init__(x0, xf, y0, yf, N, hs, i_peaks, filestr)

    def h_fun(self, x):
        # return self.h0 + np.sqrt(self.r**2 - (x-self.r)**2)
        if x <-self.r+self.dxdr or x >self.r-self.dxdr:
            return self.h0+self.r - np.sqrt(self.r**2 - (self.r-self.dxdr)**2)
        else:
            return self.h0 + self.r - np.sqrt(self.r**2 - x**2)
        
    
    