#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jan 29 08:18:07 2023

@author: sarahdennis
"""
import _graphics as graph
import numpy as np
import domain as dfd

class Height:
    
    # Takes a height function and returns discretization over given domain
    
    # h(xi) = Height.h[i]
    # h'(xi) = Height.hx[i]
    # h''(xi) = Height.hxx[i]

    def __init__(self, domain, h, h_str, h_eq, hx=None, hxx=None):
        self.h_str = h_str
        self.h_eq = h_eq
        self.hs = np.asarray([h(x) for x in domain.xs])
        
        self.h_max = max(self.hs)
        self.h_avg = np.mean(self.hs)
        
        if hx == None:
            self.hxs = dfd.center_diff(self.hs, domain)
        else:
            self.hxs = np.asarray([hx(x) for x in domain.xs])

        if hxx == None:
            self.hxxs = dfd.center_second_diff(self.hs, domain)
        else:
            self.hxxs = np.asarray([hxx(x) for x in domain.xs])
         
        
    def plot(self, domain):
        h_title = "Height (%s)"%self.h_str
        h_axis = ["Height $h(x)$", "$x$"]
        graph.plot_2D(self.hs, domain.xs, h_title, h_axis)
  
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

        super().__init__(domain, self.h, self.h_str, self.h_eq, self.hx, self.hxx)

    def h(self, x):
        return self.h0

    
    def hx(self, x):
        return 0
    
    def hxx(self, x):
        return 0

class CorrugatedHeight(Height):
    #h(x) = h_min + r(1 + cos(kx))
    def __init__(self, domain, h_mid, r, k):
        self.h_mid = h_mid
        self.r = r 
        self.k = k

        self.h_eq = "h(x) = %0.1f + %0.1f(1 + \cos(%d x))"%(self.h(domain.x0), r, k) 
        self.h_str = "Sinusoidal Height"

        super().__init__(domain, self.h, self.h_str, self.h_eq, self.hx,self.hxx)

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
        
        super().__init__(domain, self.h, self.h_str,self.h_eq, self.hx, self.hxx)


    def h(self, x):
        return self.h_min + self.m * (x - self.xf)

    def hx(self, x):
        return self.m

    def hxx(self, x):
        return 0


class StepHeight(Height):
    
    def __init__(self, domain, x_step, h_avg, r):


        self.x_step = x_step
        self.h_avg = h_avg
        self.r = r
        self.l_left = x_step-domain.x0
        self.l_right = domain.xf - x_step
        
        self.h_eq = "h(x) = {%0.1f, %0.1f}"%(h_avg + r, h_avg - r)
        self.h_str = "Rayleigh Step"
        
        super().__init__(domain, self.h, self.h_str, self.h_eq, self.hx, self.hxx)

    def h(self, x):
        if x <= self.x_step:
            return self.h_avg + self.r
        else:
            return self.h_avg - self.r

    def hx(self, x):
        return 0
    
    def hxx(self, x):
        return 0


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
        
        # For matrix solves LU, Schur, etc.
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

    # def hx(self, x):
    #     return 0
    
    # def hxx(self, x):
    #     return 0
            

    
class RandHeight(Height):
    
    def __init__(self, domain, h_max, h_min):
        self.h_max = h_max
        self.h_min = h_min
        
        self.h_str = "Random Height"
        self.h_eq = "h(x) \in [%.1f,%.1f]"%(h_min, h_max)
    
        super().__init__(domain, self.h, self.h_str, self.h_eq)

    def h(self, x):
        return (self.h_max - self.h_min) * np.random.random() + self.h_min
    

class randRectWaveHeight(Height):
    
    def __init__(self, domain, h_avg, r_max, n_steps):
        
        self.r_max = r_max
        self.h_avg = h_avg
        self.x0 = domain.x0
        self.n_steps = n_steps
        self.step_width = (domain.xf - domain.x0)/(n_steps+1)

        
        self.h_steps = np.random.uniform(low=h_avg-r_max, high=h_avg+r_max, size=n_steps+1)
        self.h_str = "%d-step Rectangle Wave"%n_steps
        self.h_eq = "h(x) \in [%.1f,%.1f]"%(h_avg-r_max, h_avg + r_max)
    
        super().__init__(domain, self.h, self.h_str, self.h_eq)
    
    def h(self, x):
        
        step_k = min(self.n_steps, int((x - self.x0)//self.step_width))
        return self.h_steps[step_k]
