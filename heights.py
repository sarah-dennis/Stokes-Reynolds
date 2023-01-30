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

    def __init__(self, domain, h, h_str, hx=None, hxx=None):
        self.h_str = h_str
        
        self.hs = [h(x) for x in domain.xs]
        
        if hx == None:
            self.hxs = dfd.center_diff(self.hs, domain)
        else:
            self.hxs = [hx(x) for x in domain.xs]
            
        if hxx == None:
            self.hxxs = dfd.center_second_diff(self.hs, domain)
    
        else:
            self.hxxs = [hxx(x) for x in domain.xs]
         
    def plot(self, domain):
        graph.plot_2D(self.hs, domain.xs, "Height %s"%self.h_str, "h" )


    def plot_all(self, domain):
        #graph.plot_2D_multi([self.hs, self.hxs, self.hxxs], domain.xs, "Height %s"%self.h_str, ["h", "hx", "hxx"])
        graph.plot_2D_multi([self.hs, self.hxs], domain.xs, "Height %s"%self.h_str, ["h", "hx"])

#------------------------------------------------------------------------------
# Example Height Functions
#------------------------------------------------------------------------------
# All example heights must have a height function h(self, x) 
# First and second derivative height functions hx(self, x) and hxx(self, x) can be specified

class CorrugatedHeight(Height):
    #h(x) = h_min + r(1 + cos(kx))
    def __init__(self, domain, h_mid, r, k):
        self.h_mid = h_mid
        self.r = r 
        self.k = k
        
        self.h_str = "h(x) = %0.1f + %0.1f(1 + \cos(%d x))"%(self.h(domain.x0), r, k) #for graph title

        super().__init__(domain, self.h, self.h_str, self.hx)

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
        self.h_max = m * (domain.x0 - domain.xf) + h_min
        
        self.h_str = "h(x) = %0.1f + %0.1f(x - %0.1f)"%(h_min, m, domain.x0) #for graph title
        
        super().__init__(domain, self.h, self.h_str, self.hx)


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
        
        
        self.h_str = "h(x) = {%0.1f, %0.1f}"%( h_left, h_right) #for graph title
        
        super().__init__(domain, self.h, self.h_str)

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
        self.h_str = "h(x) = {%0.1f, %0.1f, %0.1f}"%(h_left, h_center, h_right) #for graph title
        
        super().__init__(domain, self.h, self.h_str)

    def h(self, x):
        if x <= self.x1:
            return self.h1
        elif x <= self.x2:
            return self.h2
        else:
            return self.h3

    
    
    
    
    
    
    
    