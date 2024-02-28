 #!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jan 29 08:18:07 2023

@author: sarahdennis
"""
import graphics as graph
import numpy as np
import domain as dfd

class Height:
    
    # Takes a height function and returns discretization over given domain
    
    # h(xi) = Height.h[i]
    # h'(xi) = Height.hx[i]
    # h''(xi) = Height.hxx[i]

    def __init__(self, domain, hs, h_str, h_eq, hxs, hxxs):
        self.h_str = h_str
        self.h_eq = h_eq
        self.hs = hs
        self.hxs = hxs
        self.hxxs = hxxs
         
        self.h_max = max(self.hs)
        self.h_min = min(self.hs)
        self.h_avg = np.mean(self.hs)
        

    def plot(self, domain):
        h_title = "Height (%s)"%self.h_str
        h_axis = ["Height $h(x)$", "$x$"]
        graph.plot_2D(self.hs, domain.xs, h_title, h_axis)
  

class discreteHeight(Height):
    def __init__(self, domain, hs):
        self.h_str = "Discrete Height"
        self.h_eq = "h(x)"
        self.hs = hs
        self.hxs = dfd.center_diff(self.hs, domain)
        self.hxxs = dfd.center_second_diff(self.hs, domain)
        super().__init__(domain, self.hs, self.h_str, self.h_eq, self.hxs, self.hxxs)
#------------------------------------------------------------------------------
# Example Height Functions
#------------------------------------------------------------------------------
# All example heights must have a height function h(self, x) 
# First and second derivative height functions hx(self, x) and hxx(self, x) can be specified

class ConstantHeight(Height):
    def __init__(self, domain, h0):
        self.h_str = "Constant Height"
        self.h_eq = "h(x) = %.2f"%h0
        self.hs = np.ones(domain.Nx)*h0 
        self.hxs = np.zeros(domain.Nx)
        self.hxxs = np.zeros(domain.Nx)
        super().__init__(domain, self.hs, self.h_str, self.h_eq, self.hxs, self.hxxs)

class CorrugatedHeight(Height): #sinusoidal wave
    #h(x) = h_min + r(1 + cos(kx))
    def __init__(self, domain, h_mid, r, k):
        self.h_mid = h_mid
        self.r = r 
        self.k = k

        
        self.h_eq = "h(x) = %0.1f + %0.1f(1 + \cos(%d x))"%(self.h(domain.x0), r, k) 
        self.h_str = "Sinusoidal Height"        
        self.hs = np.asarray([self.h(x) for x in domain.xs])
        self.hxs = np.asarray([self.hx(x) for x in domain.xs])
        self.hxxs = np.asarray([self.hxx(x) for x in domain.xs])
        super().__init__(domain, self.hs, self.h_str, self.h_eq, self.hxs, self.hxxs)

    def h(self, x):
        return self.h_mid * (1 + self.r * np.cos(self.k*x))    

    def hx(self, x):
        return -self.h_mid * self.r * self.k * np.sin(self.k*x)
    
    def hxx(self, x):
        return -self.h_mid * self.r * self.k**2 * np.cos(self.k*x)


class WedgeHeight(Height): #slider bearing
    
    def __init__(self, domain, h0, h1):
        self.h0 = h0
        self.h1 = h1
        self.x0 = domain.x0
        self.m = (h1 - h0)/(domain.xf - domain.x0)
        
        self.h_eq = "h(x) = %0.1f + %0.1f(x - %0.1f)"%(self.h0, self.m, self.x0)
        self.h_str = "Wedge Slider"
        
        self.hs = np.asarray([self.h(x) for x in domain.xs])
        self.hxs = np.ones(domain.Nx)*self.m
        self.hxxs = np.zeros(domain.Nx)
        
        super().__init__(domain, self.hs, self.h_str, self.h_eq, self.hxs, self.hxxs)

    def h(self, x):
        return self.h0 + self.m * (x - self.x0)

class SawtoothHeight(Height):
    def __init__(self, domain, x_peaks, h_peaks):
        self.h_eq = "h(x)"
        self.h_str = "Piecewise Linear"
        self.h_peaks = h_peaks
        self.x_peaks = x_peaks
        self.hs, self.slopes = self.make_linear_hs(domain, x_peaks, h_peaks)
        self.hxs = dfd.center_diff(self.hs, domain)
        self.hxxs = dfd.center_second_diff(self.hs, domain)
        super().__init__(domain, self.hs, self.h_str, self.h_eq, self.hxs, self.hxxs)
        
    def make_linear_hs(self, domain, x_peaks, h_peaks):

        n_regions = len(h_peaks) - 1
        slopes = np.zeros(n_regions)

        for r in range(n_regions):
            slopes[r] = (h_peaks[r+1] - h_peaks[r])/(x_peaks[r+1] - x_peaks[r])
        
        hs = np.zeros(domain.Nx)
        r = 0
        for i in range(domain.Nx):
            xi = domain.xs[i]
            
            if xi > x_peaks[r+1] and r+1 < n_regions:

                r = r + 1
            
            
            hs[i] = h_peaks[r] + slopes[r] * (xi - x_peaks[r])
            
        return hs, slopes
    
class NStepHeight(Height): # uniform width [h1, h2, ..., hN+1]

    def __init__(self, domain, n_steps, h_steps, h_str, h_eq):
        self.n_steps = n_steps
        self.h_steps = h_steps
        self.step_width = (domain.xf - domain.x0)/(n_steps+1)

        self.hs = self.make_const_hs(domain, self.step_width, n_steps, h_steps)
        self.hxs = dfd.center_diff(self.hs, domain)
        self.hxxs = dfd.center_second_diff(self.hs, domain)
        
        super().__init__(domain, self.hs, h_str, h_eq, self.hxs, self.hxxs)

    def make_const_hs(self,domain, step_width, n_steps, h_steps):
        hs = np.zeros(domain.Nx)
        index_width = domain.Nx / (n_steps + 1)

        j=0
        for i in range(domain.Nx):

            if i >= (j+1)*index_width :
                
                j += 1
            hs[i] = h_steps[j]
 
        return hs
        
class RayleighStepHeight(NStepHeight): # N=1
    def __init__(self, domain, x_step, h1, h2):
        self.x_step = x_step

        h_str = "Step Height"
        h_eq = "h(x) = [%.2f, %.2f]"%(h1, h2)
        super().__init__(domain,1 , [h1,h2], h_str, h_eq)
    

class SquareWaveHeight(NStepHeight): #N > 1, #h(x) = h_avg +/- r
    def __init__(self, domain, h_avg, r, n_steps):
        
        h_steps = np.zeros(n_steps+1)
        
        for i in range(n_steps+1):
            h_steps[i] = h_avg + (-1)**i * r
        
        h_str = "%d-step Square Wave"%n_steps
        h_eq = "h(x) = %0.1f \pm %0.1f"%(h_avg, r)
 
        super().__init__(domain, n_steps, h_steps, h_str, h_eq)

    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    