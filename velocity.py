# -*- coding: utf-8 -*-
"""
Created on Mon Jan 26 20:21:50 2026

@author: sarah
"""

import numpy as np

import domain as dm
import reyn_velocity as rv
import reyn_velocity_ELT as rv_ELT

class Velocity:
    def __init__(self, u, v):

        self.u = u
        self.v = v
    
    def get_flux_all(self, height): #Q(x) = sum u(x) dy
        qs = np.zeros(height.Nx)
        for i in range(height.Nx):
            h = height.hs[i]
            for j in range(height.Ny):
                y = height.ys[j]    
                if y <= h:
                    qs[i] += self.u[j,i]*height.dy
                else:
                    continue
                
        return qs
    
    def get_flux(self, height):
        i=1
        Q = 0
        h = height.hs[i]
        for j in range(height.Ny):
            y = height.ys[j]    
            if y <= h:
                Q += self.u[j,i]*height.dy
            else:
                continue
        return Q
    
    def make_inc(self,height):
        u = self.u
        v = self.v
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
                        ux = dm.center_first(dx, u[j,i-1 : i+2]) #(p_E - p_W)/(2*height.dx)
                        
                    elif y < h_E and y >= h_W: # West out of bounds, fwd diff (right sided)
                        if y < h_EE:
                            ux = dm.right_first(dx, u[j, i : i+3]) #(-3*p +4*p_E -p_EE)/(2*height.dx)
                        else:
                            ux = dm.right_first_O1(dx, u[j, i: i+2]) #(p_E - p)/height.dx

                    else: # East out of bounds, bkwd diff (left sided)
                        if y < h_WW:
                            ux = dm.left_first(dx, u[j, i-2 : i+1]) #(-3*p +4*p_E -p_EE)/(2*height.dx)
                        else:
                            ux = dm.left_first_O1(dx, u[j, i-1: i+1]) 

                    if y_N < h:
                        vy  = dm.center_first(dx, v[j-1 : j+2, i])
                    else:
                        vy = dm.left_first(dx, v[j-2 : j+1, i]) #(p - p_S)/height.dx
                    
                    u_x[j,i] = ux
                    v_y[j,i] = vy

        inc = u_x + v_y
        # print('max |ux + vy|:',np.max(abs(inc)))
        return inc 
            
class ELT_Velocity(Velocity):
    
    def __init__(self,height, BC, elt_pressure):
        u, v = rv_ELT.make_ELT_velocity(height, BC, elt_pressure)
        super().__init__(u, v)
        

class Reyn_Velocity(Velocity):
    def __init__(self, height, BC, pressure):
        u, v = rv.make_reyn_velocity(height, BC, pressure)
        super().__init__(u, v)
    
        
    