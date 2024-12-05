#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec  5 13:28:42 2024

@author: sarahdennis
"""
import numpy as np

class Pressure:
    def __init__(self, height, ps):
        self.ps_1D = ps
        self.ps_2D = self.make_2D_ps(height,ps)
        
    
    def make_2D_ps(self,height,ps):
        ps_2D = np.zeros((height.Ny, height.Nx))
        
        for i in range(height.Nx):
            for j in range(height.Ny):
                
                y = height.ys[j]
                if y <= height.hs[i]:
                    ps_2D[j,i] = ps[i]
                else:
                    ps_2D[j,i] = None
                    
                
        ps_2D = np.flip(ps_2D, axis=0)
        return ps_2D
            