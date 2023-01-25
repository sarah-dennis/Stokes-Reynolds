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
            
        elif BC == 1: #fixed
            self.dx = (xf - x0)/(Nx-1)
        
        self.xs = [x0 + i*dx for i in range(Nx)]
  

class Height:
    
    def __init__(self, h, h_str, domain, hx=None):
        self.h = h
        self.grid_h = [h(x) for x in domain.xs]
        
        if hx == None:
            self.grid_hx = fd.center_diff(p)
        else:
            self.grid_hx = [hx(x) for x in domain.xs]
            
        self.string = h_str
         
    def plot(domain):
        graph.plot_2D(self.grid_height, domain.xs, "Height", "h" )









    