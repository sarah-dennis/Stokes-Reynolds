# -*- coding: utf-8 -*-
"""
Created on Wed Aug 30 12:20:23 2023

@author: sarah
"""
import numpy as np

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
                vy[j,i] = -1/(2*domain.eta) * (pressure.pxxs[i]*(1/3*domain.ys[j]**3 - 1/2*height.hs[i]*domain.ys[j]**2) - 1/2*pressure.pxs[i]*height.hxs[i]*domain.ys[j]**2) + domain.U*height.hxs[i]/(2*height.hs[i]**2)*domain.ys[j]**2
                
                #    lower surface motion:
                # vx[j,i] = 1/(2*domain.eta) * pressure.pxs[i] * (domain.ys[j]**2 - height.hs[i]*domain.ys[j]) + domain.U * (1-domain.ys[j]/height.hs[i])
                
                
                
                
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