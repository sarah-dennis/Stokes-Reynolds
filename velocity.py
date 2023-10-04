# -*- coding: utf-8 -*-
"""
Created on Wed Aug 30 12:20:23 2023

@author: sarah
"""
import numpy as np
import _graphics as graph

class Velocity:
    
    def __init__(self, vx, vy):
        
        self.vx = vx
        self.vy = vy
        
class ConstantVelocity(Velocity):
    def __init__(self, domain, height, pressure):
        domain.set_ys(height, domain.Nx)
        vx = np.zeros((domain.Ny, domain.Nx))
        vy = np.zeros((domain.Ny, domain.Nx))
        for j in range(domain.Ny):
            for i in range(domain.Nx):
                vx[j,i] = domain.U
                vy[j,i] = 0
        super().__init__(vx, vy)
 
# class CorrugatedVelocity(Velocity):
#     def __init__(self, domain, height, pressure):

class WedgeVelocity(Velocity):
    def __init__(self, domain, height, pressure):
        domain.set_ys(height, domain.Nx)
        vx = np.zeros((domain.Ny, domain.Nx))
        vy = np.zeros((domain.Ny, domain.Nx))
        for j in range(domain.Ny):
            for i in range(domain.Nx):
                #    upper surface motion:
                vx[j,i] = 1/(2*domain.eta) * pressure.pxs[i] * (domain.ys[j]**2 - height.hs[i]*domain.ys[j]) + domain.U * domain.ys[j]/height.hs[i]
    
                if i == 0 or i == domain.Nx-1:
                    vy[j,i] = 0
                else:
                    vy[j,i] = -1/(2*domain.eta) * (pressure.pxxs[i]*(1/3*domain.ys[j]**3 - 1/2*height.hs[i]*domain.ys[j]**2) - 1/2*pressure.pxs[i]*height.hxs[i]*domain.ys[j]**2) + (1/2)*domain.U*height.hxs[i]*domain.ys[j]**2/height.hs[i]**2
                
                
            
        graph.plot_2D(vy[100], domain.ys,  "Vy[10]", ["y", "vy"])
        super().__init__(vx, vy)


class SquareWaveVelocity(Velocity):
    def __init__(self, domain, height, pressure):
        domain.set_ys(height, domain.Nx)
        vx = np.zeros((domain.Ny, domain.Nx))
        vy = np.zeros((domain.Ny, domain.Nx))

        
        for j in range(domain.Ny):
            for i in range(domain.Nx):
                
                
                # upper surface motion:
                # vx[j,i] = 1/(2*domain.eta) * pressure.pxs[i] * (domain.ys[j]**2 - height.hs[i]*domain.ys[j]) + domain.U * domain.ys[j]/height.hs[i]
                
                # lower surface motion:
                vx[j,i] = 1/(2*domain.eta) * pressure.pxs[i] * (domain.ys[j]**2 - height.hs[i]*domain.ys[j]) + domain.U * (1-domain.ys[j]/height.hs[i])
                
                vy[j,i] = 0
                
        super().__init__(vx, vy)
        
        
        
class solutionVelocity(Velocity):
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
            
        super().__init__(vx, vy)
        
        
        
        
        
        