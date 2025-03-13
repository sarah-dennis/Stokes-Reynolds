#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 26 14:30:05 2025

@author: sarahdennis
"""
import numpy as np
import domain

from reyn_pressure import Pressure, FinDiffReynPressure
from reyn_velocity import Velocity, ReynoldsVelocity

class PerturbedReynSol:
    def __init__(self, height):
    

        reyn_pressure = FinDiffReynPressure(height)
        reyn_velocity = ReynoldsVelocity(height, reyn_pressure.ps_1D)
       
        
        self.x_scale = (height.xf - height.x0)/2
        self.y_scale = height.yf - height.y0
        self.Q_scale = reyn_velocity.flux
        self.P_scale = height.visc * self.Q_scale * self.x_scale * (self.y_scale**-3)
         
        self.U_scale = self.Q_scale/self.y_scale
        self.V_scale = self.Q_scale/self.x_scale                      
       
        self.u0s = reyn_velocity.vx /self.U_scale
        self.v0s = reyn_velocity.vy / self.V_scale
        self.p0s = reyn_pressure.ps_2D / self.P_scale


        delta_sqr = (self.y_scale/self.x_scale)**2
        delta_frth = (self.y_scale/self.x_scale)**4

        # make self.u2s, self.v2s, self.p2s, etc... 
        self.perturb_second(height)
        
        # make self.u4s, self.v4s, self.p4s, etc... 
        # self.perturb_fourth(height)
        
        pert_ps_2D = (self.p0s + delta_sqr * self.p2s)*self.P_scale
        pert_us_2D = (self.u0s + delta_sqr * self.u2s)*self.U_scale
        pert_vs_2D = (self.v0s + delta_sqr * self.v2s)*self.V_scale
        
        # pert_ps_2D = (self.p0s + delta_sqr * self.p2s + delta_frth * self.p4s)*self.P_scale
        # pert_us_2D = (self.u0s + delta_sqr * self.u2s + delta_frth * self.u4s)*self.U_scale
        # pert_vs_2D = (self.v0s + delta_sqr * self.v2s + delta_frth * self.v4s)*self.V_scale
        
        self.pert_pressure = Pressure(height, ps_1D = reyn_pressure.ps_1D, ps_2D=pert_ps_2D)
        self.pert_velocity = Velocity(height, pert_us_2D, pert_vs_2D)
        
        self.dP_reyn = (reyn_pressure.ps_1D[-1]-reyn_pressure.ps_1D[0])/self.P_scale
        self.dP_pert = (pert_ps_2D[0,-1]-pert_ps_2D[0,0])/self.P_scale
      
    
    def perturb_second(self, height):
        xs = height.xs/self.x_scale
        ys = height.ys/self.y_scale
        hs = height.hs/self.y_scale
        
        us = self.u0s


        h2_2xs = np.zeros(height.Nx)   # d^2/dx^2 [h^-2] @ xi
        h3_2xs = np.zeros(height.Nx)   # d^2/dx^2 [h^-3] @ xi
        h2_3xs = np.zeros(height.Nx)   # d^3/dx^3 [h^-2] @ xi
        h3_3xs = np.zeros(height.Nx)   # d^3/dx^3 [h^-3] @ xi
        
        
        # p2s = -u0_xs + c3s
        c3_xs = np.zeros(height.Nx) # d/dx [c3(x)]
        c3s = np.zeros(height.Nx) # intS d/dx[c3(x)] dx
        u0_xs = np.zeros((height.Ny, height.Nx))  # d/dx [reyn_velx] 
            
        u2s = np.zeros((height.Ny, height.Nx)) #
        v2s = np.zeros((height.Ny, height.Nx)) #

        
        dx = height.dx/(self.x_scale)
        dxx = dx**2
        dxxx = dx**3
        
        for i in range(height.Nx):
            if i >= 2  and i <= height.Nx-3 :
                h = hs[i]
                h_E = hs[i+1]                
                h_EE = hs[i+2]
                h_W = hs[i-1]
                h_WW = hs[i-2]
                
                h_x = height.hxs[i]
                
                h2_xx = ((h_E**-2) -2*(h**-2) +(h_W**-2))/dxx
                h3_xx = ((h_E**-3) -2*(h**-3) +(h_W**-3))/dxx
                
                h2_xxx = ((h_EE**-2) -2*(h_E**-2) +2*(h_W**-2) -(h_WW**-2))/(2*dxxx)
                h3_xxx = ((h_EE**-3) -2*(h_E**-3) +2*(h_W**-3) -(h_WW**-3))/(2*dxxx)
                
               
                c3_x = 6 * h2_xx * h - 18/5 * h3_xx * (h**2) 
                c3_xx = 6* (h2_xx * h_x +  h2_xxx * h) -18/5* (2*h*h_x * h3_xx + (h**2) * h3_xxx)
                
                # integrate c3_x for p2(x,y), and save {c3_x, c3} for 4th order pertubation
                c3_xs[i] = c3_x
                
                # save {h2_2x, h3_2x, h2_3x, h3_3x} for 4th order perturbation
                h2_2xs[i] = h2_xx
                h3_2xs[i] = h3_xx
                h2_3xs[i] = h2_xxx
                h3_3xs[i] = h3_xxx
                
            for j in range (height.Ny):
                y = ys[j]

                # inlet outlet buffer [x0, x1,] ... [, xN-2, xN-1]
                
                if i >= 2  and i <= height.Nx-3 and y <= h:
                    
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
                    
                    v2_A = 2 * h2_xxx * ((1/4) * (y**4) - (1/2) * (h**2) * (y**2))
                    v2_B = -2 * h2_xx * h_x * h * (y**2)
                    v2_C = -h3_xxx * ((1/5) * (y**5) - (1/2) * (h**3) * (y**2))
                    v2_D = (3/2) * h3_xx * h_x * (h**2) * (y**2)
                    v2_E = -(1/2)*c3_xx * ((1/3) * (y**3) - (1/2)* h *(y**2))
                    v2_F = (1/4) * c3_x * h_x * (y**2)
                    v2s[j,i] = v2_A + v2_B + v2_C + v2_D + v2_E + v2_F
                
                    
        # p2s[:1] = -u0_x + c3s = 0 (no correction at inlet)
        
        for i in range(2, height.Nx):
            c3s[i] =(4*c3s[i-1] -c3s[i-2] + 2*dx*c3_xs[i])/3
    
        p2s = -u0_xs + c3s  
        

        self.p2s = p2s
        self.u2s = u2s
        self.v2s = v2s
        self.c3s = c3s
        self.c3_xs = c3_xs
        
        self.h2_2xs = h2_2xs   # d^2/dx^2 [h^-2] @ xi
        self.h3_2xs = h3_2xs   # d^2/dx^2 [h^-3] @ xi
        self.h2_3xs = h2_3xs   # d^3/dx^3 [h^-2] @ xi
        self.h3_3xs = h3_3xs   # d^3/dx^3 [h^-3] @ xi
        

    def perturb_fourth(self, height):
        xs = height.xs/self.x_scale
        ys = height.ys/self.y_scale
        hs = height.hs/self.y_scale
        
        
        c5_xs = np.zeros(height.Nx) # d/dx [c5(x)]
        c5s = np.zeros(height.Nx) # intS d/dx[c5(x)] dx
        
        # p4s = -u2_xs + dxx int_0^y [v0] dy + c5s
        u2_xs = np.zeros((height.Ny, height.Nx))  # d/dx [u2] 
        v0_Sys = np.zeros((height.Ny, height.Nx)) # int_0^y [v0] dy
        v0_Sy_xxs = np.zeros((height.Ny, height.Nx)) # d^2/dx^2 int_0^y [v0] dy
        u4s = np.zeros((height.Ny, height.Nx)) #
        v4s = np.zeros((height.Ny, height.Nx)) #
        
        dy = height.dx/(self.y_scale)
        dx = height.dx/(self.x_scale)
        d2x = dx**2
        d3x = dx**3
        d4x = dx**4
        d5x = dx**5
        
        for i in range(height.Nx):
            v0_Sy_i = 0
            if i >= 3  and i <= height.Nx-4: #3-wide buffer inlet & outlet 
                
                h = hs[i]
                
                h_E = hs[i+1]                
                h_2E = hs[i+2]
                h_3E = hs[i+3]
                
                h_W = hs[i-1]
                h_2W = hs[i-2]
                h_3W = hs[i-3]

                h_x = height.hxs[i]
                
                
                h2_4x = ((h_2E**-2) -4*(h_E**-2) + 6*(h**-2) - 4*(h_W**-2) + (h_2W**-2))/d4x
                h3_4x = ((h_2E**-3) -4*(h_E**-3) + 6*(h**-3) - 4*(h_W**-3) + (h_2W**-3))/d4x
                
                h2_5x = ((h_3E**-2) -4*(h_2E**-2) + 5*(h_E**-2) - 5*(h_W**-2) + 4*(h_2W**-2) -(h_3W**-2))/(2*d5x)
                h3_5x = ((h_3E**-3) -4*(h_2E**-3) + 5*(h_E**-3) - 5*(h_W**-3) + 4*(h_2W**-3) -(h_3W**-3))/(2*d5x)
                
                # f1 = h^3 h3_xx 
                f1_2x = ((h_E**3)*self.h3_2xs[i+1] - 2*(h**3)*self.h3_2xs[i] + (h_W**3)*self.h3_2xs[i-1])/d2x
                f1_3x = ((h_2E**3)*self.h3_2xs[i+2] - 2*(h_E**3)*self.h3_2xs[i+1] + 2*(h_W**3)*self.h3_2xs[i-1] - (h_2W**3)*self.h3_2xs[i-2])/(2*d3x)
                
                # f2 = h^2 h2_xx
                f2_2x = ((h_E**2)*self.h2_2xs[i+1] - 2*(h**2)*self.h2_2xs[i] + (h_W**2)*self.h2_2xs[i-1])/d2x
                f2_3x = ((h_2E**2)*self.h2_2xs[i+2] - 2*(h_E**2)*self.h2_2xs[i+1] + 2*(h_W**2)*self.h2_2xs[i-1] - (h_2W**2)*self.h2_2xs[i-2])/(2*d3x)
                         
                # f3 = h c3_x
                f3_2x = (h_E*self.c3_xs[i+1] - 2*h*self.c3_xs[i] + h_W*self.c3_xs[i-1])/d2x
                f3_3x = (h_2E*self.c3_xs[i+2] - 2*h_E*self.c3_xs[i+1] + 2*h_W*self.c3_xs[i-1] - h_2W*self.c3_xs[i-2])/(2*d3x)

                # c3_3x = d^2/d^2x[dc3/dx]
                c3_3x = (self.c3_xs[i+1] - 2*self.c3_xs[i] + self.c3_xs[i-1])/d2x
                c3_4x = (self.c3_xs[i+2] - 2*self.c3_xs[i+1] + 2*self.c3_xs[i-1] - self.c3_xs[i-1])/(2*d3x)
                
                # c5_x = d/dx[c5]  
                c5_x = (3/14)*h3_4x*(h**4) - (3/5)*h2_4x*(h**3) - (f1_2x - 2*f2_2x)*h + (3/10)*c3_3x*(h**2) -(1/2)*f3_2x*h
                c5_2x_A = (3/14)*( h3_5x*(h**4) + 4*h3_4x*h_x*(h**3)) - (3/5)*(h2_5x*(h**3) + 3*h3_4x*h_x*(h**2))
                c5_2x_B = -((f1_3x - 2*f2_3x)*h + (f1_2x - 2*f2_2x)*h_x) + (3/10)*(c3_4x*(h**2) + 2*c3_3x*h_x*h) - (1/2)*(f3_3x*h + f3_2x*h_x) 
                c5_2x = c5_2x_A + c5_2x_B
                
                # save for p4
                c5_xs[i] = c5_x 


                for j in range (height.Ny):
                    y = ys[j]
                    v0_Sy_i += self.v0s[j,i] *  dy
                    v0_Sys[j,i] = v0_Sy_i
                    
                    # inlet outlet buffer [x0, x1,] ... [, xN-2, xN-1]
                    
                    if i >= 2  and i <= height.Nx-3 and y <= h:
                        
                        u_E = self.u2s[j,i+1]  
                        u_W = self.u2s[j,i-1]
                        x_E = xs[i+1]
                        x_W = xs[i-1]
                        
                        if y < h_E and y < h_W: # [..., i-1, i, i+1, ...]
                            u2_x = (u_E - u_W)/(2*dx)
                                    
                        elif y < h_E and y >= h_W: # [i, i+1, ...] 
            
                            x_bdry = x_E -2*dx * (h_E - y)/(h_W - h_E)
                            scale = (x_W - x_bdry)/(x_E - x_bdry)
                            u_W =  -u_E * scale
                            u2_x = (u_E - u_W)/(2*dx)
            
                        elif y >= h_E and y < h_W: # [..., i-1, i]
                            x_bdry = x_W + 2*dx * (h_W - y)/(h_E - h_W)
                            scale = (x_E - x_bdry)/(x_W -x_bdry)
                            u_E =  -u_W * scale
                            u2_x = (u_E - u_W)/(2*dx)
            
                        else: # [i]
                            u2_x = 0
                        
                        u2_xs[j,i] = u2_x

                        u4_A = (-1/20)*h3_4x*((y**6) - (h**5)*y) + (3/20)*h2_4x*((y**5) - (h**4)*y)
                        u4_B = (1/3)*(f1_2x - 2*f2_2x)*((y**3) - (h**2)*y) + (1/6)*f3_2x*((y**3) - (h**2)*y)
                        u4_C = (-1/12)*c3_3x*((y**4) - (h**3)*y) + (1/2)*c5_x*((y**2)-h*y)
                        u4s[j,i] = u4_A + u4_B + u4_C
                        
                        v4_A = (-1/20)*h3_5x*((1/7)*(y**7) - (1/2)*(h**5)*(y**2)) + (1/8)*h3_4x*h_x*(h**4)*(y**2)
                        v4_B = (3/20)*h2_5x*((1/6)*(y**6) - (1/2)*(h**4)*(y**2)) - (3/10)*h2_4x*h_x*(h**3)*(y**2)
                        v4_C = (1/3)*((f1_3x - 2*f2_3x)*((1/4)*y**4 - (1/2)*(h**2)*(y**2)) - (f1_2x - 2*f2_2x)*h_x*h*(y**2))
                        v4_D = (1/6)*(f3_3x*((1/4)*y**4 - (1/2)*(h**2)*(y**2)) - f3_2x*h_x*h*(y**2))
                        v4_E = (-1/12)*(c3_4x)*((1/5)*(y**5) - (1/2)*(h**3)*y**2) + (1/8)*c3_3x*h_x*(h**2)*(y**2)
                        v4_F = (1/2)*(c5_2x)*((1/3)*(y**3) - (1/2)*h*(y**2)) - (1/4)*c5_x*h_x*(y**2)
                        
                        v4s[j,i] = v4_A + v4_B + v4_C + v4_D + v4_E + v4_F
                        
                # p4s = -u2_xs + dxx int_0^y [v0] dy + c5s
                
                for j in range(height.Ny):
                    v0_Sy_xxs[j] = domain.center_second_diff(v0_Sys[j], height.Nx, dx)

                # p4s[:1] = 0 (no correction at inlet)
                
                for i in range(2, height.Nx):
                    c5s[i] =(4*c5s[i-1] -c5s[i-2] + 2*dx*c5_xs[i])/3
            
                p4s = -u2_xs + v0_Sy_xxs + c5s 
                
                self.p4s = p4s
                self.u4s = u4s
                self.v4s = v4s
                        
                        
                        
                        