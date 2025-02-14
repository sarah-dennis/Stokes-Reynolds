# -*- coding: utf-8 -*-
"""
Created on Wed Aug 30 12:20:23 2023

@author: sarah
"""
import numpy as np
import domain
from scipy import integrate

class Velocity:
    
    def __init__(self, height, ps):
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
                    
        for i in range(height.Nx):
            
            print(np.trapezoid(vx[:,i]), dx = height.dx)
        vy=np.flip(vy, 0)
        vx=np.flip(vx, 0)
        
        return vx, vy
    

class Adj_Velocity:
    
    def __init__(self, height, adj_ps):
        self.vx, self.vy = self.make_adj_velocity(height, adj_ps)
    

    def make_adj_velocity(self, height, ps):
        ps=np.flip(ps, 0)
        U = height.U
        visc = height.visc
        vx = np.zeros((height.Ny, height.Nx))
        vy = np.zeros((height.Ny, height.Nx))

        for j in range(height.Ny):
           
            y = height.ys[j]
            pj = ps[j]
            for i in range(1,height.Nx-1):

                h = height.hs[i]
                hx = height.hxs[i]
              
                if y <= height.hs[i]:
                    
                    px = (pj[i+1] - pj[i-1])/(2*height.dx)
                    pxx = (pj[i+1]- 2*pj[i] + pj[i-1])/(height.dx**2)
                        
                    vx[j,i] = -1/(2*visc) * px * y*(h-y) +U*(1-y/h)
                
                    a = (2*y-3*h)*pxx - 3*hx*px
                    vy[j,i] = -1/(12*visc) * (y**2) * a - (U/2)* (y/h)**2 * hx
                    
                else:
                    vx[j,i] = 0
                    vy[j,i] = 0
            
        vx[np.isnan(vx)] = 0
        flux_sum = 0
        for i in range(1,height.Nx-1):
            flux_sum += integrate.trapezoid(vx[:,i], dx=height.dx)
        flux_sum/=height.Nx-2
        print(flux_sum)
        
        vy=np.flip(vy, 0)
        vx=np.flip(vx, 0)
            
        return vx, vy
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            