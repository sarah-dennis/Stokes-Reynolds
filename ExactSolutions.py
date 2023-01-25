#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 25 18:24:34 2023

@author: sarahdennis
"""

    
class ExactSol:
    
    def __init__(self, p, domain, px=None, pxx=None):
        self.p = p
        self.grid_p = [p(x) for x in domain.xs]
        
        if px == None:
            self.grid_px = fd.center_diff(p)
        else:
            self.grid_px = [p(x) for x in domain.xs]
            
        if pxx == None:
            self.pxx = fd.center_second_diff(p)
        else:
            self.pxx = pxx
            
    def plot(domain):
        graph.plot_2D(self.grid_pressure, domain.xs, "Exact Pressure", "p" )
