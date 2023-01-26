#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 25 18:24:34 2023

@author: sarahdennis
"""

    
class Pressure:
    
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


class Corrugated(Pressure):
    
    def 
    
    def p(x):
        return p0 -6*eta*U * (h(x) + h0)/((k**2 + h0**2)*(2 + delta**2)) * h_dx(x) / h(x)**2
