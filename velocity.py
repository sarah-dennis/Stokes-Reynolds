# -*- coding: utf-8 -*-
"""
Created on Wed Aug 30 12:20:23 2023

@author: sarah
"""
import numpy as np
import _graphics as graph

class Velocity:
    
    def __init__(self, domain, height, pressure):
        self.domain = domain
        self.vel_str = "Velocity for %s (%s)"%(height.h_str, pressure.p_str)
        self.vx, self.vy = self.makeField(height, pressure)
    
    
    def makeField(self, height, pressure):
        U = self.domain.U
        eta = self.domain.eta
        self.domain.set_ys(height, self.domain.Nx)

        vx = np.zeros((self.domain.Ny, self.domain.Nx))
        vy = np.zeros((self.domain.Ny, self.domain.Nx))

        for i in range(self.domain.Nx):
            px = pressure.pxs[i]
            
            h = height.hs[i]
            hx = height.hxs[i]
            
            q = (U*h)/2 - (px*h**3)/(12*eta)

            for j in range(self.domain.Ny):
                y = self.domain.ys[j]
                
                vx[j,i] = U*(h-y)*(h-3*y)/h**2 + 6*q*y*(h-y)/h**3
                                
                vy[j,i] = -2*hx*(U/h**3 - 3*q/h**4) * y**2 * (h-y)
            
            
        mask = np.zeros((self.domain.Nx, self.domain.Ny), dtype=bool)
        for i in range(self.domain.Ny):
            for j in range(self.domain.Nx):
                mask[i,j] = self.domain.ys[i] > height.hs[j]     
        vx = np.ma.array(vx, mask=mask)
        vy = np.ma.array(vy, mask=mask)
        return vx, vy
    
      
    def plot_quivers(self):
        ax_labels = ["x", "y"]
        graph.plot_vel_quivers(self.vx, self.vy, self.domain.xs, self.domain.ys, self.vel_str, ax_labels)
        