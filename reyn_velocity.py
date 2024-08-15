# -*- coding: utf-8 -*-
"""
Created on Wed Aug 30 12:20:23 2023

@author: sarah
"""
import numpy as np
import domain

class Velocity:
    
    def __init__(self, height, ps):
        self.vel_str = "Velocity for %s "%(height.h_str)
        self.vx, self.vy = self.make_velocity(height, ps)
    
    # 2D velocity field from 1D pressure 
    def make_velocity(self, height, ps):
        U = height.U
        visc = height.visc
        vx = np.zeros((height.Ny, height.Nx))
        vy = np.zeros((height.Ny, height.Nx))

        pxs = domain.center_diff(ps, height.Nx, height.dx)
    
        for i in range(height.Nx):
            
            px = pxs[i]
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

# U = u(x,h0), Q = psi(x,h0)
# H = h0 - h(x)
def flux_to_pressure(Q, U, H):
    dp = (Q-U*H/2)*(-12/H**3)
    p0 = 0
    pN = p0 + dp
    return p0, pN

def pressure_to_flux(p0, pN, U, H):
    dp = pN - p0
    Q = dp*(H**3/-12)+U*H/2
    return Q