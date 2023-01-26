#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 25 17:33:35 2023

@author: sarahdennis
"""
import finiteDiffs as fd
import _graphics as graph

class Domain:

    def __init__(self, x0, xf, Nx, BC):
        self.x0 = x0
        self.xf = xf
        self.Nx = Nx
        
        if BC == 0: #periodic
            self.dx = (xf - x0)/(Nx)
            # note: xs[-1] = xf-dx
            
        elif BC == 1: #fixed
            self.dx = (xf - x0)/(Nx-1)
        
        self.xs = [x0 + i*dx for i in range(Nx)]
  
class Height:
    
    def __init__(self, h, h_str, domain, hx=None):
        self.h = h
        self.string = h_str

        self.grid_h = [h(x) for x in domain.xs]
        
        if hx == None:
            self.grid_hx = fd.center_diff(p)
        else:
            self.grid_hx = [hx(x) for x in domain.xs]
         
    def plot(domain):
        graph.plot_2D(self.grid_height, domain.xs, "Height", "h" )


#------------------------------------------------------------------------------
# Height Functions
#-----------------------------------------------------------------------=------
class Corrugated(Height):
     
    def __init__(self, domain, h_min, r, k):
        self.h_min = h_min
        self.r = r 
        self.k = k

        self.h_max = h_min + 2*r

        self.h_str = "h(x) = %0.1f + %0.1f(1 + \cos(%d x))"%(h0, r, k) #for graph title

        super().__init__(self.h, self.h_str, domain, self.hx)

    def h(x):
        return self.h_min * (1 + self.delta * np.cos(self.k*x))    

    def hx(x):
        return -self.h_min * self.delta * self.k * np.sin(self.k*x)


class Wedge(Height):
    
    def __init__(self, domain, h_min, m):
        self.h_min = h_min
        self.m = m
        
        self.h_max = m * (domain.x0 - domain.xf) + h_min
        
        self.h_str = "h(x) = %0.1f + %0.1f(x - %0.1f)"%(h0, m, domain.x0) #for graph title
        
        super().__init__(self.h, self.h_str, domain, self.hx)


    def h(x):
        return h_min + self.m * (x - domain.xf)

    def hx(x):
        return m




class Step(Height):
    
    def __init__(self, domain, h_left, h_right):
        self.h_left = h_left
        self.h_right = h_right
        self.x1 = (domain.xf - domain.x0)/2
        
        self.h_str = "h(x) = {%0.1f, %0.1f}"%( h_left, h_right) #for graph title
        
        super().__init__(self.h, self.h_str, domain, self.hx)

    def h(x):
        if x <= self.x1:
            return h_left
        else:
            return h_right

    
class TwoStep(Height):
    
    def __init__(self, domain, h_left, h_center, h_right):
        self.h_left = h_left
        self.h_right = h_right
        self.h_center = h_center
        self.x1 = (domain.xf - domain.x0)/3
        self.x2 = 2*(domain.xf - domain.x0)/3
        
        self.h_str = "h(x) = {%0.1f, %0.1f, %0.1f}"%(h_left, h_center, h_right) #for graph title
        
        super().__init__(self.h, self.h_str, domain, self.hx)

    def h(x):
        if x <= self.x1:
            return h_left
        elif x <= self.x2:
            return h_center
        else:
            return h_right

    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    