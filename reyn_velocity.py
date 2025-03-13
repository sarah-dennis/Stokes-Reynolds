# -*- coding: utf-8 -*-
"""
Created on Wed Aug 30 12:20:23 2023

@author: sarah
"""
import numpy as np
import domain
import graphics
from scipy import stats

class Velocity:
    def __init__(self, height, vx, vy):
        self.height=height
        self.vx = vx
        self.vy = vy
        self.flux = self.get_flux(vx)
        

            
    def get_flux(self, vx):
        lenx = vx.shape[1]
        qs = np.zeros(lenx)
    
        for i in range(self.height.Nx):
            qs[i]= np.sum(vx[:,i])*self.height.dx

        # graphics.plot_2D(qs, self.height.xs, f'flux $\lambda={self.height.lam :.2f}$', ['x','q'])
        
        q = stats.mode(qs, axis=0, keepdims=False)[0]
        return q
            
            
class ReynoldsVelocity(Velocity):
    
    def __init__(self, height, ps):
        vx, vy = self.make_velocity(height, ps)
        super().__init__(height, vx, vy)
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
                    

        # vy=np.flip(vy, 0)
        # vx=np.flip(vx, 0)
        
        return vx, vy
    

class AdjReynVelocity(Velocity):
    
    def __init__(self, height, adj_ps):
        vx, vy = self.make_adj_velocity(height, adj_ps)
        super().__init__(height, vx, vy)

    def make_adj_velocity(self, height, ps):
        # ps=np.flip(ps,0)
        vx = np.zeros((height.Ny, height.Nx))
        vy = np.zeros((height.Ny, height.Nx))

        for j in range(height.Ny):
           
            y = height.ys[j]
            pj = ps[j]
            for i in range(3, height.Nx-3):
                
                h = height.hs[i]

                p = pj[i]
                
                if y >= h:
                    vx[j,i] = 0
                    vy[j,i] = 0
                elif y == 0:
                    vx[j,i] = height.U
                    vy[j,i] = 0
                    
                else: # find px and pxx
                    h_W = height.hs[i-1]
                    p_W = pj[i-1]
                    
                    h_E = height.hs[i+1]
                    p_E = pj[i+1]
                    

                    h_WW = height.hs[i-2]
                    p_WW = pj[i-2]

                    h_WWW = height.hs[i-3]
                    p_WWW = pj[i-3]
                    

                    h_EE = height.hs[i+2]
                    p_EE = pj[i+2]


                    h_EEE = height.hs[i+3]
                    p_EEE = pj[i+3]
                    
                    if y <= h_E and y <= h_W:  # interior

                        px = (p_E - p_W)/(2*height.dx)
                        pxx = (p_E -2*p + p_W)/height.dx**2

                    elif y <= h_E and y > h_W: # West out of bounds, fwd diff (right sided)

                        if y <= h_EE:
                            px = (-3*p +4*p_E -p_EE)/(2*height.dx)
                            if y <= h_EEE:
                                pxx = (2*p -5*p_E +4*p_EE -p_EEE)/height.dx**2
                            else:
                                pxx = (p -2*p_E + p_EE)/height.dx**2

                        else:
                            px = 0#(p_E - p)/height.dx
                            pxx = 0

                    
                    elif y > h_E and y <= h_W: # East out of bounds, bkwd diff (left sided)

                        if y <= h_WW:
                            px = (3*p -4*p_W +p_WW)/(2*height.dx)
                 
                            if y <= h_WWW:
                                pxx = (2*p -5*p_W +4*p_WW -p_WWW)/height.dx**2

                            else:
                                pxx = (p -2*p_W + p_WW)/height.dx**2
                                
                        else:
                            px = (p - p_W)/height.dx
                            pxx = 0
                
                    else: # both East and West out of bounds

                        px = 0
                        pxx = 0

                    q = -height.U / h - 1/(2*height.visc) * h * px
                                        
                    qx = -1/(3*height.visc) * h * pxx

                    vx[j,i]= 1/(2*height.visc)* px * y**2 + q * y + height.U
                    
                    vy[j,i] = -1/(6*height.visc) * pxx * y**3 - 1/2 * qx * y**2
                    
                
        
        for j in range(height.Ny):
            y = height.ys[j]
            if y <= height.hs[0]:
               vx[j,0] = vx[j,3]
               vy[j,0] = vy[j,3]
            if y <= height.hs[1]:
                vx[j,1] = vx[j,3]
                vy[j,1] = vy[j,3]
            if y <= height.hs[2]:
                vx[j,2] = vx[j,3]
                vy[j,2] = vy[j,3]
            
            if y <= height.hs[-1]:
                vx[j,-1] = vx[j,height.Nx-4]
                vy[j,-1] = vy[j,height.Nx-4]
            if y <= height.hs[-2]:
                vx[j,-2] = vx[j,height.Nx-4]
                vy[j,-2] = vy[j,height.Nx-4]
            if y <= height.hs[-3]:
                vx[j,-3] = vx[j,height.Nx-4]
                vy[j,-3] = vy[j,height.Nx-4]

        # vy=np.flip(vy, 0)
        # vx=np.flip(vx, 0)
            
        return vx, vy
            
            
            
            
            