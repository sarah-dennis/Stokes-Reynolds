# -*- coding: utf-8 -*-
"""
Created on Wed Aug 30 12:20:23 2023

@author: sarah
"""
import numpy as np
import domain
import reyn_boundary as rbc
import reyn_velocity_adjusted as adj_vel
import reyn_pressure_adjusted as adj_p

class Velocity:
    def __init__(self, Q, u, v):
        # self.height=height
        self.u = u
        self.v = v
        self.Q = Q
        
    
    def get_flux(self, height):
        qs = np.zeros(height.Nx)
        for i in range(height.Nx):
            h = height.hs[i]
            for j in range(height.Ny):
                y = height.ys[j]    
                if y <= h:
                    qs[i] += self.u[j,i]*height.dy
                else:
                    continue
        for i in height.i_peaks[1:-1]:
            qs[i-3:i+4] = domain.avg_4x(qs[i-4: i+5])
            
        return qs
    
    def get_reyn_flux(self, BC, height, pressure):
        ps = pressure.ps_1D
        qs = np.zeros(height.Nx)
        
        U = BC.U
        for i in range(1,height.Nx-1):
            h = height.hs[i]
            px = domain.center_first(height.dx, ps[i-1:i+2])
            qs[i] = (U*h)/2 - (px*(h**3))/12 #/visc
        qs[0] = qs[1]
        qs[-1] = qs[-2]
        return qs
    
    def get_adj_flux(self, BC, height, pressure):
        reyn_ps = pressure.ps_1D
        qs = np.zeros(height.Nx)
        U = BC.U
        
        
        pxs = domain.center_diff(reyn_ps, height.Nx, height.dx)
        p2xs = domain.center_second_diff(reyn_ps, height.Nx, height.dx)
        p3xs = domain.center_third_diff(reyn_ps, height.Nx, height.dx)
        p4xs = domain.center_fourth_diff(reyn_ps, height.Nx, height.dx)
        
        for i in height.i_peaks[1:-1]:
            p2xs[i-1:i+2] =domain.avg_2x(p2xs[i-2 : i+3]) 
            p3xs[i-2:i+3] =domain.avg_3x(p3xs[i-3 : i+4]) 
            p4xs[i-2:i+3] = domain.avg_3x(p4xs[i-3 : i+4])
        sigmas, sxs, s2xs = adj_p.make_sigmas(height,BC, pxs,p2xs,p3xs,p4xs)
        
        for i in range(height.Nx):
            

            h = height.hs[i]
            hx = height.hxs[i]
            h2x = height.h2xs[i]
            px = pxs[i]
            p2x = p2xs[i]
            p3x = p3xs[i]
            sx = sxs[i]
            
            q_re = (U*h)/2 - (px*(h**3))/12 #/visc
            q_padj = p3x*(h**5)/80 -(h*p3x +2*p2x*hx + px*h2x)*(h**4)/48
            q_uadj = U*(h2x*(h**2)-2*(hx**2)*h)/24 
            qs[i] = q_re + q_padj + q_uadj - sx*(h**3)/12
            #/visc  
        for i in height.i_peaks[1:-1]:
            qs[i-2:i+3] = domain.avg_3x(qs[i-3 : i+4])
            
        return qs
    
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
    def __init__(self, height, BC, ps=None) :
        if isinstance(BC, rbc.Fixed):
            h0=height.hs[0]
            px0 = domain.center_first(height.dx, ps[0:3])
            Q = (BC.U*h0)/2 - (px0*(h0**3))/12 #/visc
            
        elif isinstance(BC, rbc.Mixed):
            Q = BC.Q
            
        u, v = self.make_velocity(height, BC.U, Q)
        super().__init__(Q, u, v)
        
    # 2D velocity field from 1D pressure 
    def make_velocity(self, height, U, Q):
        # visc = height.visc
        u = np.zeros((height.Ny, height.Nx))
        v = np.zeros((height.Ny, height.Nx))


        for i in range(height.Nx):

            h = height.hs[i]
            hx = height.hxs[i]

            for j in range(height.Ny):
                y = height.ys[j]
                if y <= height.hs[i]:
                    u[j,i] = (h-y)*(U*(h-3*y)/h**2 + 6*Q*y/h**3)
                    v[j,i] = -2*hx * y**2 * (h-y) *(U/h**3 - 3*Q/h**4)
                    
                else:
                    u[j,i] = 0
                    v[j,i] = 0
                    
        return u, v

class TGAdj_ReynVelocity(Velocity):
    
    def __init__(self, height, BC, adj_pressure):
        u, v = adj_vel.make_adj_velocity_TG(height, BC, adj_pressure)
        if isinstance(BC, rbc.Mixed):
            Q = BC.Q
        else:
            h0=height.hs[0]
            px0 = domain.right_first(height.dx, adj_pressure.ps_2D[0,0:3])
            Q = (BC.U*h0)/2 - (px0*(h0**3))/12 #/visc
        super().__init__(Q, u, v)
     
class VelAdj_ReynVelocity(Velocity):
    
    def __init__(self, height, BC, adj_pressure):
        u, v = adj_vel.make_adj_velocity(height, BC, adj_pressure)
        if isinstance(BC, rbc.Mixed):
            Q = BC.Q
        else:
            h0=height.hs[0]
            px0 = domain.right_first(height.dx, adj_pressure.ps_2D[0,0:3])
            Q = (BC.U*h0)/2 - (px0*(h0**3))/12 #/visc
        super().__init__(Q, u, v)
        
