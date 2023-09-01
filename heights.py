#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jan 29 08:18:07 2023

@author: sarahdennis
"""
import _graphics as graph
import numpy as np
import domain as dfd
from matplotlib import pyplot as pp

class Height:
    
    # Takes a height function and returns discretization over given domain
    
    # h(xi) = Height.h[i]
    # h'(xi) = Height.hx[i]
    # h''(xi) = Height.hxx[i]

    def __init__(self, domain, h, h_str, h_eq, hx=None, hxx=None):
        self.h_str = h_str
        self.h_eq = h_eq
        
        self.hs = [h(x) for x in domain.xs]
        
        self.h_max = max(self.hs)
        self.ys = np.linspace(0, 1.1* self.h_max, domain.Nx)
        
        if hx == None:
            self.hxs = dfd.center_diff(self.hs, domain)
            
        else:
            
            self.hxs = np.asarray([hx(x) for x in domain.xs])

        if hxx == None:
            self.hxxs = dfd.center_second_diff(self.hs, domain)
    
        else:
            self.hxxs = np.asarray([hxx(x) for x in domain.xs])
         
    def plot(self, domain):
        fig = pp.figure()
        pp.plot(domain.xs, self.hs, color='b')

        pp.title(self.h_str)
        
        pp.xlabel("$x$")
        pp.ylabel("Height $h(x)$")
        
        pp.ylim(0, 1.2*max(self.hs))
        pp.xticks(np.arange(domain.xs[0], domain.xs[-1]+0.001, self.step_width))
     
        return fig
        


    def plot_derivs(self, domain):
        graph.plot_2D_multi([self.hs, self.hxs, self.hxxs], domain.xs, "Height %s"%self.h_str, ["h", "hx", "hxx"])
        #graph.plot_2D_multi([self.hs, self.hxs], domain.xs, "Height %s"%self.h_str, ["h", "hx"])

#------------------------------------------------------------------------------
# Example Height Functions
#------------------------------------------------------------------------------
# All example heights must have a height function h(self, x) 
# First and second derivative height functions hx(self, x) and hxx(self, x) can be specified
class ConstantHeight(Height):
    def __init__(self, domain, h0):
        self.h0 = h0
        self.h_str = "Constant Height"
        self.h_eq = "h(x) = %.2f"%h0
        super().__init__(domain, self.h, self.h_str, self.h_eq)
        
    def h(self, x):
        return self.h0
    
class CorrugatedHeight(Height):
    #h(x) = h_min + r(1 + cos(kx))
    def __init__(self, domain, h_mid, r, k):
        self.h_mid = h_mid
        self.r = r 
        self.k = k
        
        self.h_eq = "h(x) = %0.1f + %0.1f(1 + \cos(%d x))"%(self.h(domain.x0), r, k) 
        self.h_str = "Sinuosoidal Height"

        super().__init__(domain, self.h, self.h_str, self.h_eq, self.hx)

    def h(self, x):
        return self.h_mid * (1 + self.r * np.cos(self.k*x))    

    def hx(self, x):
        return -self.h_mid * self.r * self.k * np.sin(self.k*x)
    
    def hxx(self, x):
        return -self.h_mid * self.r * self.k**2 * np.cos(self.k*x)


class WedgeHeight(Height):
    
    def __init__(self, domain, h_min, m):
        self.h_min = h_min
        self.m = m
        self.xf = domain.xf
        
        self.h_eq = "h(x) = %0.1f + %0.1f(x - %0.1f)"%(h_min, m, domain.x0)
        self.h_str = "Wedge Slider"
        
        super().__init__(domain, self.h, self.h_str,self.h_eq, self.hx)


    def h(self, x):
        return self.h_min + self.m * (x - self.xf)

    def hx(self, x):
        return self.m

    def hxx(self, x):
        return 0

class StepHeight(Height):
    
    def __init__(self, domain, h_left, h_right):
        self.h_left = h_left
        self.h_right = h_right
        self.x1 = (domain.xf - domain.x0)/2
        
        self.l_left = self.x1-domain.x0
        self.l_right = domain.xf - self.x1
        
        self.h_eq = "h(x) = {%0.1f, %0.1f}"%( h_left, h_right)
        self.h_str = "Rayleigh Step"
        
        super().__init__(domain, self.h, self.h_str, self.h_eq)

    def h(self, x):
        if x <= self.x1:
            return self.h_left
        else:
            return self.h_right

class TwoStepHeight(Height):
    
    def __init__(self, domain, h_left, h_center, h_right):
        self.h1 = h_left
        self.h2 = h_center
        self.h3 = h_right
        self.x1 = (domain.xf - domain.x0)/3
        self.x2 = 2*(domain.xf - domain.x0)/3
        self.l1 = self.x1 - domain.x0
        self.l2 = self.x2 - self.x1
        self.l3 = domain.xf - self.x2
        self.h_str = "Rayleigh Two Step"
        self.h_eq = "h(x) = {%0.1f, %0.1f, %0.1f}"%(h_left, h_center, h_right) 
        
        super().__init__(domain, self.h, self.h_str, self.h_eq)

    def h(self, x):
        if x <= self.x1:
            return self.h1
        elif x <= self.x2:
            return self.h2
        else:
            return self.h3


class SquareWaveHeight(Height):
    
    def __init__(self, domain, h_avg, r, n_steps):
        
        self.r = r
        self.h_avg = h_avg
        self.n_steps = n_steps
        self.step_width = (domain.xf - domain.x0)/(n_steps+1)
        
        if 0.01 * (h_avg + r) < r:
            # (Li & Chen, 2007) : roughness height < 0.01 * total height 
            print("Warning: r > 0.01 * (h + r)")
            print("%.5f > %.5f"%(r, 0.01 * (h_avg + r)))
        
        self.h_steps = np.zeros(n_steps+1)
       
        for i in range(n_steps+1):
            x = domain.x0 + self.step_width * (i + 0.5)
            self.h_steps[i] = self.h(x)
        
        self.h_str = "%d-step Square Wave"%n_steps
        self.h_eq = "h(x) = %0.1f \pm %0.1f"%(h_avg, r)
 
        super().__init__(domain, self.h, self.h_str, self.h_eq)
        
    def h(self, x):
        if np.sin(np.pi * x/self.step_width) >= 0:
            return self.h_avg + self.r
        else:
            return self.h_avg - self.r
        
            
        
    
        

    
    
    
    
    
    
    