# -*- coding: utf-8 -*-
"""
Created on Wed Aug 30 12:20:23 2023

@author: sarah
"""
import numpy as np
from numpy.polynomial import Polynomial as poly
import _graphics as graph

class Velocity:
    
    def __init__(self, vx, vy):
        
        self.vx = vx
        self.vy = vy
        
    def plot_vx_y0(self, domain, j):
        vx_title = "Velocity $Vx$ at $y_0=%.2f$"%domain.ys[j]
        vx_labels = ["$x$", "Velocity $Vx(x, y_0)$"]
        graph.plot_2D(self.vx[j], domain.xs, vx_title, vx_labels)
        c = poly.fit(domain.xs, self.vx[j], deg = 2)
        print("Vx(x, %.1f): %s"%(domain.ys[j],c))
        
        
    def plot_vx_x0(self, domain, i):
        vx_title = "Velocity $Vx$ at $x_0=%.2f$"%domain.xs[i]
        vx_labels = ["Velocity $Vx(x_0, y)$","$y$"]
        graph.plot_2D(domain.ys, self.vx[:,i], vx_title, vx_labels)

        c = poly.fit(domain.ys, self.vx[:,i], deg = 2)
        print("Vx(%.1f, y): %s"%(domain.xs[i],c))
        
        
        
    def plot_vy_y0(self, domain, j):
        vy_title = "Velocity $Vy$ at $x_0=%.2f$"%domain.ys[j]
        vy_labels = ["$x$", "Velocity $Vy(x, y_0)$"]
        graph.plot_2D(self.vy[j], domain.xs, vy_title, vy_labels)
    
    def plot_vy_x0(self, domain, i):
        vy_title = "Velocity $Vy$ at $x_0=%.2f$"%domain.xs[i]
        vy_labels = ["Velocity $Vy(x_0, y)$","$y$" ]
        graph.plot_2D(domain.ys, self.vy[:,i], vy_title, vy_labels) 



class SolutionVelocity(Velocity):
    def __init__(self, domain, height, pressure):
        U = domain.U
        eta = domain.eta
        domain.set_ys(height, domain.Nx)

        vx = np.zeros((domain.Ny, domain.Nx))
        vy = np.zeros((domain.Ny, domain.Nx))

        for i in range(domain.Nx):
            px = pressure.pxs[i]
            
            h = height.hs[i]
            hx = height.hxs[i]
            
            q = (U*h)/2 - (px*h**3)/(12*eta)

            for j in range(domain.Ny):
                y = domain.ys[j]
                
                vx[j,i] = U*(h-y)*(h-3*y)/h**2 + 6*q*y*(h-y)/h**3
                                
                vy[j,i] = -2*hx*(U/h**3 - 3*q/h**4) * y**2 * (h-y)
            
            
        mask = np.zeros((domain.Nx, domain.Ny), dtype=bool)
        for i in range(domain.Ny):
            for j in range(domain.Nx):
                mask[i,j] = domain.ys[i] > height.hs[j]     
        vx = np.ma.array(vx, mask=mask)
        vy = np.ma.array(vy, mask=mask)        
        
        super().__init__(vx, vy)
        
        
        
        
        
        