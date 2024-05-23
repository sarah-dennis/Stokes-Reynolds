# -*- coding: utf-8 -*-
"""
Created on Tue May 21 15:15:15 2024

@author: sarah
"""
import numpy as np
from domain import Height

class triangle(Height):
    def __init__(self, x0, xf, y0, yf, N, U, Re, filestr):
        # N even --> (xL-x0)/2 is triangle vertex
        # slope dividing 2N  --> maximal true boundary points
        # self.x0 = x0
        # self.xf = xf
        # self.y0 = y0
        # self.yf = yf
        self.slope = int(2*yf/xf)
        self.filename = filestr
        # self.N = N 
        # self.h = 1/N
        # self.n = (xf-x0)*N + 1
        # self.m = (yf-y0)*N + 1
        
        Nx = (xf-x0)*N + 1
        hs = np.zeros(Nx) # Reynolds needs this, Stokes its built into dPsi
        # TODO: make hs since it could be helpful to make for plotting
        
        super().__init__(x0, xf, y0, yf, N, hs, U, Re, filestr)