#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 26 14:30:05 2025

@author: sarahdennis
"""
import numpy as np
import domain
import graphics

from reyn_pressure import Pressure, FinDiffReynPressure
from reyn_velocity import Velocity, ReynoldsVelocity

class PerturbedReynSol:
    def __init__(self, height, order):
        self.order = order
        if order < 0 or order > 4:
            return Exception(f"order {order} not in range [0,4]")
            
        self.reyn_pressure = FinDiffReynPressure(height)
        self.reyn_velocity = ReynoldsVelocity(height, self.reyn_pressure.ps_1D)
    
        
        
        self.x_scale = (height.xf - height.x0)/2
        self.y_scale = height.yf - height.y0
        self.Q_scale = self.reyn_velocity.flux

        self.P_scale = height.visc * self.Q_scale * self.x_scale * (self.y_scale**-3)
         
        self.U_scale = self.Q_scale/self.y_scale
        self.V_scale = self.Q_scale/self.x_scale                      

        self.u0s = self.reyn_velocity.vx /self.U_scale
        self.v0s = self.reyn_velocity.vy / self.V_scale
        self.p0s = self.reyn_pressure.ps_2D / self.P_scale
        
        self.dP_reyn = (self.p0s[0,-1]-self.p0s[0,0])
        
        delta = self.y_scale/self.x_scale
        if order > 1:

            # make self.u2s, self.v2s, self.p2s, etc... 
            self.perturb_second(height)
            pert2_ps_2D = (self.p0s + (delta**2) * self.p2s)*self.P_scale
            pert2_us_2D = (self.u0s + (delta**2) * self.u2s)*self.U_scale
            pert2_vs_2D = (self.v0s + (delta**2) * self.v2s)*self.V_scale
            self.pert2_pressure = Pressure(height, ps_1D = self.reyn_pressure.ps_1D, ps_2D=pert2_ps_2D)
            self.pert2_velocity = Velocity(height, pert2_us_2D, pert2_vs_2D)
            self.dP_pert2 = (self.p2s[0,-1]-self.p2s[0,0])

        if order > 2: 

            # make self.u4s, self.v4s, self.p4s, etc... 
            self.perturb_fourth(height)

            pert4_ps_2D = pert2_ps_2D + (delta**4) * self.p4s*self.P_scale
            pert4_us_2D = pert2_us_2D + (delta**4) * self.u4s*self.U_scale
            pert4_vs_2D = pert2_vs_2D + (delta**4) * self.v4s*self.V_scale
        
            self.pert4_pressure = Pressure(height, ps_1D = self.reyn_pressure.ps_1D, ps_2D=pert4_ps_2D)
            self.pert4_velocity = Velocity(height, pert4_us_2D, pert4_vs_2D)
        
            self.dP_pert4 = (self.p4s[0,-1]-self.p4s[0,0])
    
    
        
    def perturb_second(self, height):

        ys = height.ys/self.y_scale
        hs = height.hs/self.y_scale
        
        h_xs = np.zeros(height.Nx)     # d/dx[h] @ xi
        h2_2xs = np.zeros(height.Nx)   # d^2/dx^2 [h^-2] @ xi
        h3_2xs = np.zeros(height.Nx)   # d^2/dx^2 [h^-3] @ xi
        h2_3xs = np.zeros(height.Nx)   # d^3/dx^3 [h^-2] @ xi
        h3_3xs = np.zeros(height.Nx)   # d^3/dx^3 [h^-3] @ xi
        
        
        # p2s = -u0_xs + c3s
        c3_2xs = np.zeros(height.Nx)                # d^2/dx^2 [c3(x)]
        c3_xs = np.zeros(height.Nx)                 # d/dx [c3(x)]
        c3s = np.zeros(height.Nx)                   # intS d/dx[c3(x)] dx
        u0_xs = np.zeros((height.Ny, height.Nx))    # d/dx [reyn_velx] 
            
        u2s = np.zeros((height.Ny, height.Nx)) #
        v2s = np.zeros((height.Ny, height.Nx)) #

        
        dx = height.dx/(self.x_scale)
        dxx = dx**2
        dxxx = dx**3
        
        for i in range(height.Nx):
            h = hs[i]
            if i < 2:
                h_E = hs[i+1] 
                h_EE = hs[i+2] 
                
                if i == 0:
                    h_W = h_E
                    h_WW = h_EE
    
                elif i == 1:
                    h_W = hs[i-1]
                    h_WW = h

            elif i > height.Nx-3:
                h_W = hs[i-1]
                h_WW = hs[i-2]
                
                if i == height.Nx-1:
                    h_E = h_W
                    h_EE = h_WW
                    
                elif i == height.Nx-2:
                    h_E = hs[i+1]
                    h_EE = h
                
            else:  
                
                h_E = hs[i+1] 
                h_EE = hs[i+2]
                h_W = hs[i-1]
                h_WW = hs[i-2]
                
            h_x = (h_E - h_W)/(2*dx)
            
            h2_xx = ((h_E**-2) -2*(h**-2) +(h_W**-2))/dxx
            h3_xx = ((h_E**-3) -2*(h**-3) +(h_W**-3))/dxx
            
            h2_xxx = ((h_EE**-2) -2*(h_E**-2) +2*(h_W**-2) -(h_WW**-2))/(2*dxxx)
            h3_xxx = ((h_EE**-3) -2*(h_E**-3) +2*(h_W**-3) -(h_WW**-3))/(2*dxxx)
            
            c3_x  = 6 * h2_xx * h - 18/5 * h3_xx * (h**2) 
            c3_xx = 6 * (h2_xx * h_x +  h2_xxx * h) - 18/5* (2*h*h_x * h3_xx + (h**2) * h3_xxx)

            # integrate c3_x for p2(x,y), and save {c3_x, c3} for 4th order pertubation
            c3_xs[i]  = c3_x
            c3_2xs[i] = c3_xx
            
            # save {h2_2x, h3_2x, h2_3x, h3_3x} for 4th order perturbation
            h_xs[i]   = h_x
            h2_2xs[i] = h2_xx
            h3_2xs[i] = h3_xx
            h2_3xs[i] = h2_xxx
            h3_3xs[i] = h3_xxx
            
            for j in range (height.Ny):
                y = ys[j]
    
                if y <= h:
                    u0_xs[j,i] = 6*h_x * (-2*y*(h**-3) + 3*(y**2)* (h**-4))

                    u2_A = -2 * h2_xx * (y**3 - (h**2) * y) 
                    u2_B = h3_xx * (y**4 - (h**3) * y)
                    u2_C = (1/2) * c3_x * (y**2 - h * y)
                    u2s[j,i] = (u2_A + u2_B +  u2_C)
                    
                    v2_A = h2_xxx * ((1/2) * (y**4) - (h**2) * (y**2))
                    v2_B = -2 * h2_xx * h_x * h * (y**2)
                    v2_C = -h3_xxx * ((1/5) * (y**5) - (1/2) * (h**3) * (y**2))
                    v2_D = (3/2) * h3_xx * h_x * (h**2) * (y**2)
                    v2_E = -(1/2)*c3_xx * ((1/3) * (y**3) - (1/2)* h *(y**2))
                    v2_F = (1/4) * c3_x * h_x * (y**2)
                    v2s[j,i] = (v2_A + v2_B + v2_C + v2_D + v2_E + v2_F)
            

                
                    
        # p2s[:1] = -u0_x + c3s = 0 (no correction at inlet)
        
        c3s[0] = np.mean(u0_xs[:,0]) #u0_xs[0]
        c3s[1] = c3s[0] + c3_xs[2]*dx
        for i in range(2, height.Nx):
            c3s[i] =(4*c3s[i-1] -c3s[i-2] + 2*dx*c3_xs[i])/3
    
        p2s = -u0_xs + c3s 

        self.p2s = p2s
        self.u2s = u2s
        self.v2s = v2s
        
        
        self.c3s = c3s
        self.c3_xs = c3_xs
        self.c3_2xs = c3_2xs
        self.h_xs = h_xs
        self.h2_2xs = h2_2xs   # d^2/dx^2 [h^-2] @ xi
        self.h3_2xs = h3_2xs   # d^2/dx^2 [h^-3] @ xi
        self.h2_3xs = h2_3xs   # d^3/dx^3 [h^-2] @ xi
        self.h3_3xs = h3_3xs   # d^3/dx^3 [h^-3] @ xi
        

    def perturb_fourth(self, height):

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
        
        dy = height.dy/(self.y_scale)
        dx = height.dx/(self.x_scale)
        d2x = dx**2
        d3x = dx**3
        
        for i in range(height.Nx):
            v0_Sy_i = 0
            
            h = hs[i]
            
            h_x = self.h_xs[i]
            h2_2x = self.h2_2xs[i]
            h3_2x = self.h3_2xs[i]
            h2_3x = self.h2_3xs[i]
            h3_3x = self.h3_3xs[i]
            
            c3_x = self.c3_xs[i]
            c3_2x = self.c3_2xs[i]
            
            if i < 2: # Inlet  boundary 
            
                h_E = hs[i+1] 
                h_2E = hs[i+2] 

                h2_2x_E = self.h2_2xs[i+1]
                h2_2x_2E = self.h2_2xs[i+2]
                h3_2x_E = self.h3_2xs[i+1]
                h3_2x_2E = self.h3_2xs[i+2]
                h2_3x_E = self.h2_3xs[i+1]
                h3_3x_E = self.h3_3xs[i+1]
                
                c3_x_E = self.c3_xs[i+1]
                c3_x_2E = self.c3_xs[i+2]
                c3_2x_E = self.c3_2xs[i+1]
                
                if i == 0: # adjust W:=i-1 & 2W:=i-2
                    h_W = h_E
                    h_2W = h_2E
                    
                    h2_2x_W = h2_2x_E
                    h2_2x_2W = h2_2x_2E
                    h3_2x_W = h3_2x_E
                    h3_2x_2W = h3_2x_2E
                    h2_3x_W = h2_3x_E
                    h3_3x_W = h3_3x_E
                    
                    c3_x_W = c3_x_E
                    c3_x_2W = c3_x_2E
                    c3_2x_W = c3_2x_E
                    
    
                elif i == 1: # adjust 2W:=i-2
                    h_W = hs[i-1]
                    h_2W = h
                    
                    h2_2x_W = self.h2_2xs[i-1]
                    h2_2x_2W = h2_2x
                    h3_2x_W = self.h3_2xs[i-1]
                    h3_2x_2W = h3_2x
                    h2_3x_W = self.h2_3xs[i-1]
                    h3_3x_W = self.h3_3xs[i-1]
                    
                    c3_x_W = self.c3_xs[i-1]
                    c3_x_2W = c3_x
                    c3_2x_W = self.c3_2xs[i-1]
                    

            if i > height.Nx-3: # Outlet boundary 
                h_W = hs[i-1]
                h_2W = hs[i-2]

                h2_2x_W = self.h2_2xs[i-1]
                h2_2x_2W = self.h2_2xs[i-2]
                h3_2x_W = self.h3_2xs[i-1]
                h3_2x_2W = self.h3_2xs[i-2]
                h2_3x_W = self.h2_3xs[i-1]
                h3_3x_W = self.h3_3xs[i-1]
                
                c3_x_W = self.c3_xs[i-1]
                c3_x_2W = self.c3_xs[i-2]
                c3_2x_W = self.c3_2xs[i-1]
                
                if i == height.Nx-1: # adjust E:=i+1 & 2E:=i+2
                    h_E = h_W
                    h_2E = h_2W
                
                    h2_2x_E = h2_2x_W
                    h2_2x_2E = h2_2x_2W
                    h3_2x_E = h3_2x_W
                    h3_2x_2E = h3_2x_2W
                    h2_3x_E = h2_3x_W
                    h3_3x_E = h3_3x_W
                    
                    c3_x_E = c3_x_W
                    c3_x_2E = c3_x_2W
                    c3_2x_E = c3_2x_W
                    
                elif i == height.Nx-2: # adjust 2E:=i+2
                    h_E = hs[i+1]
                    h_2E = h_E
                    
                    h2_2x_E = self.h2_2xs[i+1]
                    h2_2x_2E = h2_2x
                    h3_2x_E = self.h3_2xs[i+1]
                    h3_2x_2E = h3_2x
                    h2_3x_E = self.h2_3xs[i+1]
                    h3_3x_E = self.h3_3xs[i+1]
                    
                    
                    c3_x_E = self.c3_xs[i+1]
                    c3_x_2E = c3_x
                    c3_2x_E = self.c3_2xs[i+1]
                
            else: # all interior nbrs
                
                h_E  = hs[i+1] 
                h_2E = hs[i+2]
                h_W  = hs[i-1]
                h_2W = hs[i-2]
                
                h2_2x_E  = self.h2_2xs[i+1]
                h2_2x_2E = self.h2_2xs[i+2]
                h3_2x_E  = self.h3_2xs[i+1]
                h3_2x_2E = self.h3_2xs[i+2]
                h2_3x_E  = self.h2_3xs[i+1]
                h3_3x_E  = self.h3_3xs[i+1]
                
                h2_2x_W  = self.h2_2xs[i-1]
                h2_2x_2W = self.h2_2xs[i-2]
                h3_2x_W  = self.h3_2xs[i-1]
                h3_2x_2W = self.h3_2xs[i-2]
                h2_3x_W  = self.h2_3xs[i-1]
                h3_3x_W  = self.h3_3xs[i-1]
                
                c3_x_E  = self.c3_xs[i+1]
                c3_x_2E = self.c3_xs[i+2]
                c3_2x_E = self.c3_2xs[i+1]
                
                c3_x_W  = self.c3_xs[i-1]
                c3_x_2W = self.c3_xs[i-2]
                c3_2x_W = self.c3_2xs[i-1]
        
            
        

            h2_4x = (h2_3x_E - h2_3x_W)/(2*dx)
            h2_5x = (h2_3x_E - 2*h2_3x + h2_3x_W)/d2x
            
            h3_4x = (h3_3x_E - h3_3x_W)/(2*dx)
            h3_5x = (h3_3x_E - 2*h3_3x + h3_3x_W)/d2x
            
            
            f1    = (h**3)   * h3_2x
            f1_E  = (h_E**3) * h3_2x_E    
            f1_2E = (h_2E**3)* h3_2x_2E
            f1_W  = (h_W**3) * h3_2x_W
            f1_2W = (h_2W**3)* h3_2x_2W                
            
            f1_2x = (f1_E - 2*f1 + f1_W)/d2x
            f1_3x = (f1_2E - 2*f1_E + 2*f1_W - f1_2W)/(2*d3x)
            
            f2    = (h**2)   * h2_2x
            f2_E  = (h_E**2) * h2_2x_E
            f2_2E = (h_2E**2)* h2_2x_2E
            
            f2_W  = (h_W**2) * h2_2x_W
            f2_2W = (h_2W**2)* h2_2x_2W
            
            f2_2x = (f2_E - 2*f2 + f2_W)/d2x
            f2_3x = (f2_2E - 2*f2_E + 2*f2_W - f2_2W)/(2*d3x)
            
            f3    = h    * c3_x
            f3_E  = h_E  * c3_x_E
            f3_2E = h_2E * c3_x_2E
            f3_W  = h_W  * c3_x_W
            f3_2W = h_2W * c3_x_2W #fixed typo 3/18
            
            f3_2x = (f3_E - 2*f3 + f3_W)/d2x
            f3_3x = (f3_2E - 2*f3_E + 2*f3_W - f3_2W)/(2*d3x)


            c3_3x = (c3_2x_E - c3_2x_W)/(2*dx)
            c3_4x = (c3_2x_E - 2*c3_2x + c3_2x_W)/d2x
            
            c5_x = (3/14)*h3_4x*(h**4) - (3/5)*h2_4x*(h**3) - (f1_2x - 2*f2_2x)*h + (3/10)*c3_3x*(h**2) -(1/2)*f3_2x*h
            
            c5_2x_A = (3/14)*(h3_5x*(h**4) + 4*h3_4x*h_x*(h**3)) 
            c5_2x_B = -(3/5)*(h2_5x*(h**3) + 3*h2_4x*h_x*(h**2))
            c5_2x_C = -((f1_3x - 2*f2_3x)*h + (f1_2x - 2*f2_2x)*h_x) 
            c5_2x_D = (3/10)*(c3_4x*(h**2) + 2*c3_3x*h_x*h) 
            c5_2x_E = -(1/2)*(f3_3x*h + f3_2x*h_x) 
            c5_2x = c5_2x_A + c5_2x_B + c5_2x_C + c5_2x_D + c5_2x_E 


            if i <3 or i > height.Nx-4:
                print(h2_4x, h2_5x, h3_4x, h3_5x, c3_3x, c3_4x, c5_x, c5_2x)    
            # save for p4
            c5_xs[i] = c5_x

            for j in range (height.Ny):
                y = ys[j]
                
                v0_Sy_i += self.v0s[j,i] *  dy
                
                v0_Sys[j,i] = v0_Sy_i
                
                if y <= h:
         
                    u2x_A = -2*h2_3x*((y**3)-(h**2)*y) + 4*h2_2x*h_x*h*y
                    u2x_B = h3_3x*((y**4)-(h**3)*y) -3*h3_2x*h_x*(h**2)*y
                    u2x_C = (1/2)*(c3_2x *((y**2)-h*y) - h_x*y) 
                    u2_xs[j,i] = u2x_A + u2x_B + u2x_C
                    
                    u4_A = (-1/20)*h3_4x*((y**6) - (h**5)*y) + (3/20)*h2_4x*((y**5) - (h**4)*y)
                    u4_B = (1/3)*(f1_2x - 2*f2_2x)*((y**3) - (h**2)*y) + (1/6)*f3_2x*((y**3) - (h**2)*y)
                    u4_C = (-1/12)*c3_3x*((y**4) - (h**3)*y) + (1/2)*c5_x*((y**2)-h*y)
                    u4s[j,i] = (u4_A + u4_B + u4_C)
                    
                    v4_A = (-1/20)*h3_5x*((1/7)*(y**7) - (1/2)*(h**5)*(y**2)) + (1/8)*h3_4x*h_x*(h**4)*(y**2)
                    v4_B = (3/20)*h2_5x*((1/6)*(y**6) - (1/2)*(h**4)*(y**2)) - (3/10)*h2_4x*h_x*(h**3)*(y**2)
                    v4_C = (1/3)*((f1_3x - 2*f2_3x)*((1/4)*(y**4) - (1/2)*(h**2)*(y**2)) - (f1_2x - 2*f2_2x)*h_x*h*(y**2))
                    v4_D = (1/6)*(f3_3x*((1/4)*(y**4) - (1/2)*(h**2)*(y**2)) - f3_2x*h_x*h*(y**2))
                    v4_E = (-1/12)*(c3_4x)*((1/5)*(y**5) - (1/2)*(h**3)*(y**2)) + (1/8)*c3_3x*h_x*(h**2)*(y**2)
                    v4_F = (1/2)*(c5_2x)*((1/3)*(y**3) - (1/2)*h*(y**2)) - (1/4)*c5_x*h_x*(y**2)
                    v4s[j,i] = -(v4_A + v4_B + v4_C + v4_D + v4_E + v4_F)
    

                        
        # p4s = -u2_xs + dxx int_0^y [v0] dy + c5s
                
        for j in range(height.Ny):
            v0_Sy_xxs[j] = domain.center_second_diff(v0_Sys[j], height.Nx, dx)

        p4s = np.zeros((height.Ny, height.Nx))
        
        c5s[0] = 0
        c5s[1] = c5s[0] + c5_xs[1]*dx

        for i in range(2, height.Nx):
            c5s[i] =(4*c5s[i-1] -c5s[i-2] + 2*dx*c5_xs[i])/3

        p4s = -u2_xs + v0_Sy_xxs + c5s 

        # p4s[:,-5:]=0
        self.p4s = p4s
        self.u4s = u4s
        self.v4s = v4s
                        
                        
                        
                        