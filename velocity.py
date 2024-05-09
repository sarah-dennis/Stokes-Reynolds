# -*- coding: utf-8 -*-
"""
Created on Wed Aug 30 12:20:23 2023

@author: sarah
"""
import numpy as np

class Velocity:
    
    def __init__(self, height, pressure):
        self.vel_str = "Velocity for %s (%s)"%(height.h_str, pressure.p_str)
        self.vx, self.vy = self.make_velocity(height, pressure)
    
    # 2D velocity field from 1D pressure under Reynolds assumptions
    def make_velocity(self, height, pressure):
        U = height.U
        visc = height.visc

        vx = np.zeros((height.Ny, height.Nx))
        vy = np.zeros((height.Ny, height.Nx))

        for i in range(height.Nx):
            px = pressure.pxs[i]
            
            h = height.hs[i]
            hx = height.hxs[i]
            
            q = (U*h)/2 - (px*h**3)/(12*visc)

            for j in range(height.Ny):
                y = height.ys[j]
                if y <= height.hs[i]:
                    vx[j,i] = U*(h-y)*(h-3*y)/h**2 + 6*q*y*(h-y)/h**3
                                
                    vy[j,i] = -2*hx*(U/h**3 - 3*q/h**4) * y**2 * (h-y)
                else:
                    vx[j,i] = 0
                                
                    vy[j,i] = 0
        

        return vx, vy

