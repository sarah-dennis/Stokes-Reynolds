# -*- coding: utf-8 -*-
"""
Created on Wed Aug 30 12:20:23 2023

@author: sarah
"""
import numpy as np
import domain
import graphics
from scipy import stats
import reyn_velocity_adjusted as adj_vel

class Velocity:
    def __init__(self, height, vx, vy):
        self.height=height
        self.vx = vx
        self.vy = vy
        self.flux, self.qs = self.get_flux(vx)

        self.inc = self.make_inc(height, vx, vy)

    def get_flux(self, vx):
        lenx = vx.shape[1]
        qs = np.zeros(lenx)
    
        for i in range(self.height.Nx):
            h = self.height.hs[i]
            for j in range (self.height.Ny):
                
                    y = self.height.ys[j]    
                
                    if y <= h:
                        qs[i]+= vx[j,i]*self.height.dy
                    else:
                        continue

        q = qs[0]
        return q, qs
    
    
    def make_inc(self, height, u, v):
        inc = np.zeros((height.Ny, height.Nx))
        u_x = np.zeros((height.Ny, height.Nx))
        v_y = np.zeros((height.Ny, height.Nx))
        dx = height.dx
        for j in range(3, height.Ny-3):
           
            y = height.ys[j]
            
            for i in range(2, height.Nx-2):
                
                h = height.hs[i]
                
                if y >= h:
                    u_x[j,i] = 0
                    v_y[j,i] = 0
                
                else: # find u_x and v_y
                    h_W = height.hs[i-1]
                    h_E = height.hs[i+1]
                    h_WW = height.hs[i-2]
                    h_EE = height.hs[i+2]
                    y_N = height.ys[j+1]                 
                    if y < h_E and y < h_W:  # interior
                        ux = domain.center_first(dx, u[j,i-1 : i+2]) #(p_E - p_W)/(2*height.dx)
                        
                    elif y < h_E and y >= h_W: # West out of bounds, fwd diff (right sided)
                        if y < h_EE:
                            ux = domain.right_first(dx, u[j, i : i+3]) #(-3*p +4*p_E -p_EE)/(2*height.dx)
                        else:
                            ux = domain.right_first_O1(dx, u[j, i: i+2]) #(p_E - p)/height.dx

                    else: # East out of bounds, bkwd diff (left sided)
                        if y < h_WW:
                            ux = domain.left_first(dx, u[j, i-2 : i+1]) #(-3*p +4*p_E -p_EE)/(2*height.dx)
                        else:
                            ux = domain.left_first_O1(dx, u[j, i-1: i+1]) 

                    if y_N < h:
                        vy  = domain.center_first(dx, v[j-1 : j+2, i])
                    else:
                        vy = domain.left_first(dx, v[j-2 : j+1, i]) #(p - p_S)/height.dx
                    
                    u_x[j,i] = ux
                    v_y[j,i] = vy

        inc = u_x + v_y
        # print('max |ux + vy|:',np.max(abs(inc)))
        return inc 
            
class ReynVelocity(Velocity):
    def __init__(self, height, ps):
        vx, vy = self.make_velocity(height, ps)
        super().__init__(height, vx, vy)
        
    # 2D velocity field from 1D pressure 
    def make_velocity(self, height, ps):
        U = height.U
        visc = height.visc
        vx = np.zeros((height.Ny, height.Nx))
        vy = np.zeros((height.Ny, height.Nx))

        h0=height.hs[0]
        px0 = domain.right_first(height.dx, ps[0:3])
        q = (U*h0)/2 - (px0*(h0**3))/(12*visc)
        # print('Reyn flux:', q)

        for i in range(height.Nx):

            h = height.hs[i]
            hx = height.hxs[i]

            for j in range(height.Ny):
                y = height.ys[j]
                if y <= height.hs[i]:
                    vx[j,i] = (h-y)*(U*(h-3*y)/h**2 + 6*q*y/h**3)
                    vy[j,i] = -2*hx * y**2 * (h-y) *(U/h**3 - 3*q/h**4)
                    
                else:
                    vx[j,i] = 0
                    vy[j,i] = 0
                    
        return vx, vy

class Adjusted_ReynVelocity_TG(Velocity):
    
    def __init__(self, height, adj_pressure):
        vx, vy = adj_vel.make_adj_velocity_TG(height, adj_pressure)
        super().__init__(height, vx, vy)
        
# class Adjusted_ReynVelocity_inc(Velocity):
    
#     def __init__(self, height, adj_ps):
#         vx, vy = adj_vel.make_adj_velocity(height, adj_ps)
#         super().__init__(height, vx, vy)
        
class Adjusted_ReynVelocity(Velocity):
    
    def __init__(self, height, adj_pressure):
        vx, vy = adj_vel.make_adj_velocity(height, adj_pressure)
        super().__init__(height, vx, vy)
        
