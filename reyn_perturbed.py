#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 26 14:30:05 2025

@author: sarahdennis
"""
import numpy as np

from reyn_pressure import Pressure, FinDiffReynPressure
from reyn_velocity import Velocity, ReynoldsVelocity

class PerturbedReynSol:
    def __init__(self, height):
    

        self.reyn_pressure = FinDiffReynPressure(height)
        self.reyn_velocity = ReynoldsVelocity(height, self.reyn_pressure.ps_1D)
        
        self.x_scale = (height.xf - height.x0)/2
        self.y_scale = height.yf - height.y0
        self.Q_scale = self.reyn_velocity.flux

        self.P_scale = height.visc * self.Q_scale * self.x_scale * (self.y_scale**-3)
         
        self.U_scale = self.Q_scale/self.y_scale
        self.V_scale = self.Q_scale/self.x_scale                      
       

        delta_sqr = (self.y_scale/self.x_scale)**2

        
        p2s, u2s, v2s = self.perturb_second(height)
        

        
        pert_ps_2D = self.reyn_pressure.ps_2D + delta_sqr * p2s*self.P_scale
        pert_us_2D = self.reyn_velocity.vx + delta_sqr * u2s*self.U_scale
        pert_vs_2D = self.reyn_velocity.vy + delta_sqr * v2s*self.V_scale
        
        self.pert_pressure = Pressure(height, ps_1D = self.reyn_pressure.ps_1D, ps_2D=pert_ps_2D)
        self.pert_velocity = Velocity(height, pert_us_2D, pert_vs_2D)
        
        self.dP_reyn = (self.reyn_pressure.ps_1D[-1]-self.reyn_pressure.ps_1D[0])/self.P_scale
        self.dP_pert = (pert_ps_2D[0,-1]-pert_ps_2D[0,0])/self.P_scale
      
    
    def perturb_second(self, height):
        xs = height.xs/self.x_scale
        ys = height.ys/self.y_scale
        hs = height.hs/self.y_scale
        us = self.reyn_velocity.vx/self.U_scale

        # p2s = -u0_xs + c3s
        c3_xs = np.zeros(height.Nx) # d/dx [c3(x)]
        c3s = np.zeros(height.Nx) # intS d/dx[c3(x)] dx
        u0_xs = np.zeros((height.Ny, height.Nx))  # d/dx [reyn_velx] 
            
        u2s = np.zeros((height.Ny, height.Nx)) #
        v2s = np.zeros((height.Ny, height.Nx)) #
        
        # h_x = 0     # d/dx [h] @ xi
        # h2_xx = 0   # d^2/dx^2 [h^-2] @ xi
        # h3_xx = 0   # d^2/dx^2 [h^-3] @ xi
        # h2_xxx = 0  # d^3/dx^3 [h^-2] @ xi
        # h3_xxx = 0  # d^3/dx^3 [h^-3] @ xi
        
        dx = height.dx/(self.x_scale)
        dxx = dx**2
        dxxx = dx**3
        
        for i in range(height.Nx):
            if i < 2  or i > height.Nx-3:
                c3_xs[i] = 0
                
            else:
                h = hs[i]
                h_E = hs[i+1]                
                h_EE = hs[i+2]
                h_W = hs[i-1]
                h_WW = hs[i-2]
                
                h_x = (h_E - h_W)/(2*dx)    
                
                
                h2_xx = ((h_E**-2) -2*(h**-2) +(h_W**-2))/dxx
                h3_xx = ((h_E**-3) -2*(h**-3) +(h_W**-3))/dxx
                 
                h2_xxx = ((h_EE**-2) -2*(h_E**-2) +2*(h_W**-2) -(h_WW**-2))/(2*dxxx)
                h3_xxx = ((h_EE**-3) -2*(h_E**-3) +2*(h_W**-3) -(h_WW**-3))/(2*dxxx)
               
                # print(h2_xx, h2_xx_b)
                
                c3_x = 6 * h2_xx * h - 18/5 * h3_xx * (h**2) 
                c3_xx = 6* (h2_xx * h_x +  h2_xxx * h) -18/5* (2*h*h_x * h3_xx + (h**2) * h3_xxx)
                
                c3_xs[i] = c3_x
                
            for j in range (height.Ny):
                y = ys[j]

                # inlet outlet buffer [x0, x1,] ... [, xN-2, xN-1]
                if i < 2  or i > height.Nx-3 or y >= h:
                    u0_xs[j,i] = 0
                    u2s[j,i] = 0
                    v2s[j,i] = 0
                
                else:
                    u_E = us[j,i+1]  
                    u_W = us[j,i-1]
                    x_E = xs[i+1]
                    x_W = xs[i-1]
                    
                    if y < h_E and y < h_W: # [..., i-1, i, i+1, ...]
                        u0_x = (u_E - u_W)/(2*dx)
                                
                    elif y < h_E and y >= h_W: # [i, i+1, ...] 

                        x_bdry = x_E -2*dx * (h_E - y)/(h_W - h_E)
                        scale = (x_W - x_bdry)/(x_E - x_bdry)
                        u_W =  -u_E * scale
                        u0_x = (u_E - u_W)/(2*dx)

                        
                    elif y >= h_E and y < h_W: # [..., i-1, i]
                        x_bdry = x_W + 2*dx * (h_W - y)/(h_E - h_W)
                        scale = (x_E - x_bdry)/(x_W -x_bdry)
                        u_E =  -u_W * scale
                        u0_x = (u_E - u_W)/(2*dx)

                           
                    else: # [i]
                        u0_x = 0
                    
                    u0_xs[j,i] = u0_x
                    
                    u2_A = -2 * h2_xx * (y**3 - (h**2) * y) 
                    u2_B = h3_xx * (y**4 - (h**3) * y)
                    u2_C = (1/2) * c3_x * (y**2 - h * y)
                    u2s[j,i] = u2_A + u2_B +  u2_C
                    
                    v2_A = h2_xxx * ((1/2) * (y**4) - (h**2) * (y**2))
                    v2_B = -2 * h2_xx * h_x * (y**2)
                    v2_C = -h3_xxx * ((1/5)*(y**5) - (h**3) * (y**2))
                    v2_D = (3/2) * h3_xx * (h_x**2) * y**2
                    v2_E = -c3_xx * ((1/6)*(y**3) - (1/4)* h *(y**2))
                    v2_F = (1/4) * c3_x * h_x * (y**2)
                    v2s[j,i] = v2_A + v2_B + v2_C + v2_D + v2_E + v2_F
                
                    
        c3s[0] = 0 # p2s[0] = 0 (no correction at inlet)
        c3s[1] = 0
        for i in range(2, height.Nx):
            c3s[i] =(4*c3s[i-1] -c3s[i-2] + 2*dx*c3_xs[i])/3
    
        p2s = -u0_xs + c3s  
        

    
        return p2s, u2s, v2s









