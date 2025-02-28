#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 26 14:30:05 2025

@author: sarahdennis
"""
import numpy as np
from scipy import integrate



def get_pert2_p_u_v(height, reyn_ps, reyn_velx):
    x_scale, y_scale, U_scale, V_scale, P_scale = height.get_dimless_vars()
    # P2(x,y) = -dU0/dx + C3(x)
    
    # xs = height.xs/x_scale
    ys = height.ys/y_scale
    hs = height.hs/y_scale
    reyn_ps /= P_scale


    u0_xs = np.zeros((height.Ny, height.Nx))  # d/dx [reyn_velx] 
    c3x_xs = np.zeros((height.Ny, height.Nx)) # d/dx [dC3/dx]
    
    u2s = np.zeros((height.Ny, height.Nx)) #
    v2s = np.zeros((height.Ny, height.Nx)) #
    
    h_x=0       # d/dx [h] @ xi
    h2_xx = 0   # d^2/dx^2 [h^-2] @ xi
    h3_xx = 0   # d^2/dx^2 [h^-3] @ xi
    h2_xxx = 0  # d^3/dx^3 [h^-2] @ xi
    h3_xxx = 0  # d^3/dx^3 [h^-3] @ xi
    
    h2_xx_E = 0 # d^2/dx^2 [h^-2] @ i+1
    h2_xx_W = 0 # d^2/dx^2 [h^-2] @ i-1
    h3_xx_E = 0 # d^2/dx^2 [h^-3] @ i+1
    h3_xx_W = 0 # d^2/dx^2 [h^-3] @ i-1
    
    dx = height.dx/x_scale
    dy = height.dx/y_scale
    dxx = dx**2
    dxxx = dx**3
    
    
    for j in range (height.Ny):
        y = ys[j]
        
        for i in range(height.Nx):
            h = hs[i] 
            
            # finite difference buffer [x0, x1, x2] ... [xN-3, xN-2, xN-1]
            if i < 3  or i > height.Nx-4 or y >= h:
                u0_xs[j,i] = 0
                c3x_xs[j,i] = 0
                u2s[j,i] = 0
                v2s[j,i] = 0
            
            else:
                
                h_E = height.hs[i+1]                
                h_EE = height.hs[i+2]
                h_EEE = height.hs[i+3]
                h_W = height.hs[i-1]
                h_WW = height.hs[i-2]
                h_WWW = height.hs[i-3]
                
                u = reyn_velx[j,i]
                u_E = reyn_velx[j,i+1]
                u_EE = reyn_velx[j,i+2]   
                u_W = reyn_velx[j,i-1]
                u_WW = reyn_velx[j,i+2]

                
                if y < h_E and y < h_W: # [..., i-1, i, i+1, ...]
                    
                    u0_x = (u_E - u_W)/(2*dx)
                    
                    h_x = (h_E - h_W)/(2*dx)
                    
                    h2_xx = ((h_E**-2) -2*(h**-2) + (h_W**-2))/(dxx)
                    h3_xx = ((h_E**-3) -2*(h**-3) + (h_W**-3))/(dxx)
            
                    
                    if y < h_EE and y < h_WW: # [i-2, i-1, i, i+1, i+2]
                        
                        h2_xxx = ((h_EE**-2) - 2*(h_E**-2) +2*(h_W**-2) - (h_WW**-2))/(2*dxxx)
                        h3_xxx = ((h_EE**-3) - 2*(h_E**-3) +2*(h_W**-3) - (h_WW**-3))/(2*dxxx)
                        
                        h2_xx_E = ((h_EE**-2) -2*(h_E**-2) + (h**-2))/dxx
                        h2_xx_W = ((h**-2)    -2*(h_W**-2) + (h_WW**-2))/dxx

                        h3_xx_E = ((h_EE**-3) -2*(h_E**-3) + (h**-3))/dxx
                        h3_xx_W = ((h**-3)    -2*(h_W**-3) + (h_WW**-3))/dxx                  
                    
                    else: # one or both of i+2 or i-2 is out of bounds
                        
                        h2_xxx = 0
                        h3_xxx = 0
                        
                        if y < h_EE and y>= h_WW: # [i-1, i, i+1, i+2]
                       
                            h2_xx_E = ((h_EE**-2) -2*(h_E**-2) + (h**-2))/dxx
                            h2_xx_W = ((h_E**-2)  -2*(h**-2)   + (h_W**-2))/dxx
                            
                            h3_xx_E = ((h_EE**-3) -2*(h_E**-3) + (h**-3))/dxx
                            h3_xx_W = ((h_E**-3)  -2*(h**-3)   + (h_W**-3))/dxx                        
                            
                        elif y >= h_EE and y < h_WW: #  [i-2, i-1, i, i+1]
                            
                            h2_xx_E = ((h_E**-2) -2*(h**-2)   + (h_W**-2))/dxx
                            h2_xx_W = ((h**-2)   -2*(h_W**-2) + (h_WW**-2))/dxx
          
                            h3_xx_E = ((h_E**-3) -2*(h**-3)   + (h_W**-3))/dxx
                            h2_xx_W = ((h**-3)   -2*(h_W**-3) + (h_WW**-3))/dxx

                        else: # [i-1, i, i+1]
                            h2_xx_E = 0
                            h2_xx_W = 0
                            h3_xx_E = 0
                            h3_xx_W = 0
                    
                    
                    c3x_E = 6 * h_E * h2_xx_E - 18/5 * (h_E**2) * h3_xx_E
                    c3x_W = 6 * h_W * h2_xx_W - 18/5 * (h_W**2) * h3_xx_W
                    
                    c3x_x = (c3x_E - c3x_W)(2*height.dx)
                    
                            
                elif y < h_E and y >= h_W: # [i, i+1, ...] 
                    
                    h2_xxx = 0
                    h3_xxx = 0
                                        
                    h2_xx_W = 0
                    h3_xx_W = 0
                    
                    if y < h_EE: # [i, i+1, i+2, ...]
                    
                        u0_x = (-3*u +4*u_E -u_EE)/(2*dx)   
                        h_x = (h_EE - h)/(2*dx)
                        h2_xx_E = ((h_EE**-2) -2*(h_E**-2) + (h**-2))/dxx
                        h3_xx_E = ((h_EE**-3) -2*(h_E**-3) + (h**-3))/dxx

                        if y < h_EEE: #[i, i+1, i+2, i+3]
                            h2_xx = (2*(h**-2) -5*(h_E**-2) +4*(h_EE**-2) -(h_EEE**-2))/dxx
                            h3_xx = (2*(h**-3) -5*(h_E**-3) +4*(h_EE**-3) -(h_EEE**-3))/dxx
                            
                        else: # [i, i+1, i+2]
                            h2_xx = ((h_EE**-2) - 2*(h_E**-2) + (h**-2))/dxx
                            h3_xx = ((h_EE**-3) - 2*(h_E**-3) + (h**-3))/dxx
            
                    else: # [i, i+1]
                        u0_x = (u_E - u)/dx
                        h_x = (h_E - h)/dx
                        h2_xx = 0
                        h3_xx = 0
                        h2_xx_E = 0
                        h3_xx_E = 0
                    
                    c3x_E = 6 * h_E * h2_xx_E - 18/5 * (h_E**2) * h3_xx_E
                    c3x   = 6 * h   * h2_xx   - 18/5 * (h**2)   * h3_xx
                    
                    c3x_x = (c3x_E - c3x)/dx

                elif y >= h_E and y < h_W: # [..., i-1, i]
                    h2_xxx = 0
                    h3_xxx = 0
                                        
                    h2_xx_E = 0
                    h3_xx_E = 0

                    if y < h_WW: # [..., i-2, i-1, i]
                        u0_x = (3*u -4*u_W +u_WW)/(2*dx)
                        h_x = (h - h_WW)/(2*dx)
                        h2_xx_W = ((h_WW**-2) -2*(h_W**-2) + (h**-2))/dxx
                        h3_xx_W = ((h_WW**-3) -2*(h_W**-3) + (h**-3))/dxx
                        
                        if y < h_WWW: # [i-3, i-2, i-1, i]
                            h2_xx = (2*(h**-2) -5*(h_W**-2) +4*(h_WW**-2) -(h_WWW**-2))/dxx
                            h3_xx = (2*(h**-3) -5*(h_W**-3) +4*(h_WW**-3) -(h_WWW**-3))/dxx 
                            
                        else: # [i-2, i-1, i]
                            h2_xx = ((h**-2) - 2*(h_W**-2) + (h_WW**-2))/dxx
                            h3_xx = ((h**-3) - 2*(h_W**-3) + (h_WW**-3))/dxx
            
                    else: # [i-1, i]
                        u0_x = (u - u_W)/dx
                        h2_xx = 0
                        h3_xx = 0
                        h2_xx_W = 0
                        h3_xx_W = 0
                    
                    c3x_W = 6 * h_W * h2_xx_W - 18/5 * (h_W**2) * h3_xx_W
                    c3x = 6 * h * h2_xx - 18/5 * (h**2) * h3_xx
                    
                    c3x_x = (c3x - c3x_W)/dx

                else: # [i]
                    u0_x = 0
                    h2_xx = 0
                    h3_xx = 0
                    h2_xxx = 0
                    h3_xxx = 0
                    c3x_x = 0
                
                u0_xs[j,i] = u0_x
                c3x_xs[j,i] = c3x_x
                
                u2s[j,i] = -2 * h2_xx * (y**3 - (h**2) * y) + h3_xx * (y**4 - (h**3) * y) + (1/2) * c3x_xs[j,i] * (y**2 - h * y)      
                v2s[j,i] = h2_xxx * (0.5 * (y**4) - (h**2) * (y**2)) - 2 * h2_xx * h_x * (y**2) - h3_xxx * ((1/5)*y**5 - (h**3) * (y**2)) + (3/2) * h3_xx * (h_x**2) * y**2 - c3x_xs[j,i]*((1/6)*(y**3)-(1/4)*h*(y**2)) + (1/4)*c3x * h_x * (y**2)
    
    c3s = np.zeros(height.Nx)
    c3_init = u0_xs[height.Ny//2,0]
    for i in range(height.Nx):
        c3s[i] = integrate.trapezoid(c3x_xs[:,i], x=ys, dx=dy) + c3_init
    
    
    p2s = -u0_xs + c3s    
    
    p2s_reDim = p2s*P_scale
    u2s_reDim = u2s*U_scale
    v2s_reDim = v2s*V_scale


    return p2s_reDim, u2s_reDim, v2s_reDim










