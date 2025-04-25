# -*- coding: utf-8 -*-
"""
Created on Wed Aug 30 12:20:23 2023

@author: sarah
"""
import numpy as np
import domain
import graphics

class Velocity:
    def __init__(self, height, vx, vy):
        self.height=height
        self.vx = vx
        self.vy = vy
        self.flux = self.get_flux(vx)
        
        
        self.inc = self.make_inc(height, vx, vy)

    def get_flux(self, vx):
        lenx = vx.shape[1]
        qs = np.zeros(lenx)
    
        for i in range(2, self.height.Nx-2):
            h = self.height.hs[i]
            for j in range (self.height.Ny):
                
                    y = self.height.ys[j]    
                
                    if y <= h:
                        qs[i]+= vx[j,i]*self.height.dx
                    else:
                        continue

        graphics.plot_2D(qs[3:-3], self.height.xs[3:-3], 'flux', ['x', 'q'])
        q = np.average(qs)
        # print(q)
        return q
    
    
    
    
    def make_inc(self, height, u, v):
        inc = np.zeros((height.Ny, height.Nx))
        u_x = np.zeros((height.Ny, height.Nx))
        v_y = np.zeros((height.Ny, height.Nx))
        dx = height.dx
        for j in range(3, height.Ny-3):
           
            y = height.ys[j]
            
            for i in range(4, height.Nx-4):
                
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
        # print(np.max(abs(inc)))
        return inc 
            
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

        h0=height.hs[0]
        px0 = domain.right_first(height.dx, ps[0:3])
        q = (U*h0)/2 - (px0*(h0**3))/(12*visc)

        for i in range(height.Nx):

            h = height.hs[i]
            hx = height.hxs[i]

            for j in range(height.Ny):
                y = height.ys[j]
                if y <= height.hs[i]:
                    vx[j,i] = U*(h-y)*(h-3*y)/h**2 + 6*q*y*(h-y)/h**3
                    vy[j,i] = -2*hx*(U/h**3 - 3*q/h**4) * y**2 * (h-y)
                    
                else:
                    vx[j,i] = 0
                    vy[j,i] = 0
                    
        return vx, vy
    

class AdjReynVelocity(Velocity):
    
    def __init__(self, height, adj_ps):
        vx, vy = self.make_adj_velocity(height, adj_ps)
        
        super().__init__(height, vx, vy)
        
    def make_adj_velocity(self, height, ps):

        
        # Find v(x,y) = int_y0^y[p] dy - int_y0^h [p] dy * (y-y0)/h
        vs = self.make_vs(height, ps)

        # u(x,y) = int_x0^x[du/dx] dx
        us = self.make_us(height, vs, ps)
       

        return us, vs
    
    def make_intPdys(self, height, ps):
        intPdys = np.zeros((height.Ny, height.Nx)) # int_y0^y[p] dy
        for i in range(height.Nx):
            h = height.hs[i]
            intPdy_i = 0
            
            for j in range(1,height.Ny):
                y = height.ys[j]
        
                if y <= h:
                    intPdy_i += ps[j-1,i] *height.dy
                    intPdys[j,i] = intPdy_i
                else:
                    continue
        # graphics.plot_contour_mesh(intPdys, height.xs, height.ys, '$\int p dy$', ['int p dy','x','y'], vmin=-1, vmax=1)
        return intPdys
    
    def make_vs(self, height, ps):
        vs = np.zeros((height.Ny, height.Nx))
        intPdys = self.make_intPdys(height, ps)
        
        for i in range(height.Nx):
            h = height.hs[i]
            jh = int((h/height.dy))-1
            for j in range(height.Ny):
                y = height.ys[j]
                if y <= h:
                    vs[j,i] =  intPdys[j,i] - intPdys[jh,i]*y/h
                else:
                    continue
        return vs
                    
    def make_dvdys(self, height, vs):
        vys = np.zeros((height.Ny,height.Nx))
        # dv/dy = -du/dx
        for i in range(height.Nx):
            h = height.hs[i]
            # [-vs[1,i], vs[0,i], vs[1,i]]
            vys[0,i]= domain.right_first(height.dy,vs[0:3,i])
            # vyB = domain.center_first(height.dy, [-vs[1,i], vs[0,i], vs[1,i]])
            
            for j in range(1,height.Ny):
                y = height.ys[j]    
                if y < h:
                    if j < height.Ny-1 and height.ys[j+1] <= h:
                        vys[j,i] = domain.center_first(height.dy,vs[j-1:j+2,i])
                        
                    elif j < height.Ny-2 and height.ys[j+2] <= h:
                        vys[j,i] = domain.left_first(height.dy,vs[j-2:j+1,i])
                    # else: 
                        # vys[j,i] = 0  
                elif y == h: 
                    vys[j,i] = domain.left_first(height.dy,vs[j-2:j+1,i])
                    continue
                # else:
                    # vys[j,i] = 0  
                    
        graphics.plot_contour_mesh(vys, height.xs, height.ys, 'vys', ['vys','x','y'], vmin=-1, vmax=1)

        return vys
    
    def make_us(self, height, vs, ps):
        vys = self.make_dvdys(height, vs)

        h0 = height.hs[0]
        px0 = domain.right_first(height.dx, ps[0,0:3])
        q = (height.U*h0)/2 - (px0*(h0**3))/(12*height.visc)
        
        us = np.zeros((height.Ny, height.Nx))
        
        for i in range(height.Nx):
            h = height.hs[i]
            x = height.xs[i]
            
            for j in range(height.Ny):
                
                y = height.ys[j]
                if j == 0:
                    us[j,i] = height.U
                    
                elif y == h:
                    us[j,i] = 0 
                    continue # y > h

                elif y < h:
                
                    if i == 0:
                        us[j,i] = height.U*(h-y)*(h-3*y)/h**2 + 6*q*y*(h-y)/h**3 
        
                    else:
         
                        if y <= height.hs[i-1]: 
                            
                            if i > 1 and y <= height.hs[i-2]: # contour using i-1 and i-2
                                us[j,i] = (4*us[j,i-1] -us[j,i-2] -2*vys[j,i]*height.dx)/3
                                
                            else: # contour using i-1 only
                                us[j,i] = us[j,i-1] - vys[j,i] * height.dx
                            
                        else: # (i-1,j) oob
                            
                            x_bdry_W = height.xs[i+1] + height.dx * (y - height.hs[i+1])/(h - height.hs[i+1])
                            us[j,i] = - vys[j,i] * (x-x_bdry_W)
        return us
    
    # def make_adj_velocity_old(self, height, ps):

    #     dx = height.dx
    #     hs = height.hs
    #     hxs = height.hxs
    #     ys = height.ys
    #     us = np.zeros((height.Ny, height.Nx))
    #     vs = np.zeros((height.Ny, height.Nx))
    #     pxs = np.zeros((height.Ny, height.Nx))
    #     pxxs = np.zeros((height.Ny, height.Nx))
    #     px_hs = np.zeros(height.Nx)
    #     pxx_hs = np.zeros(height.Nx)
    #     f1_x_err = np.zeros(height.Nx)
    
    #     for j in range(height.Ny):
           
    #         y = ys[j]
    #         pj = ps[j]
            
    #         for i in range(3,height.Nx-3):
                 
               
    #             h = hs[i]

    #             if y <= h:
    #                 if y <= hs[i+1] and y <= hs[i-1]: #interior nbrs
    #                     px = domain.center_first(dx, pj[i-1 : i+2])
    #                     pxx = domain.center_second(dx, pj[i-1 : i+2]) 
                    
    #                 elif y <= hs[i+1] and y > hs[i-1]: # West out of bounds, fwd diff (right sided)
    #                     if y <= hs[i+2]:
    #                         px = domain.right_first(dx, pj[i : i+3]) 
                            
    #                         if y <= hs[i+3]:
    #                             pxx = domain.right_second(dx, pj[i : i+4])
                                
    #                         else:
    #                             pxx = domain.right_second_O1(dx, pj[i : i+3]) 
    
    #                     else:
    #                         px = domain.right_first_O1(dx, pj[i: i+2]) 
    #                         pxx = 0
                    
    #                 elif y > hs[i+1] and y <= hs[i-1]: # East out of bounds, bkwd diff (left sided)
    #                     if y <= hs[i-2]:
    #                         px = domain.left_first(dx, pj[i-2 : i+1]) 
                            
    #                         if y <= hs[i-3]:
    #                             pxx = domain.left_second(dx, pj[i-3 : i+1]) 
    #                         else:
    #                             pxx = domain.left_second_O1(dx, pj[i-2 : i+1])
                                
    #                     else:
    #                         px = domain.left_first_O1(dx, pj[i-1 : i+1]) 
    #                         pxx = 0
                           
    #                 else: # both East and West out of bounds
    #                     px = 0
    #                     pxx = 0
                
    #             else: #exterior
    #                 px = 0
    #                 pxx = 0
                
    #             pxs[j,i] = px
    #             pxxs[j,i] = pxx
                
    #             if y == h or (y < h and y+dx >= h):
    #                 px_hs[i] = px
    #                 pxx_hs[i] = pxx

    #     for i in range(height.Nx):
    #         h = hs[i]
    #         hx = hxs[i]
            
            
    #         f1 =  -1/(2*height.visc) * h * px_hs[i] - height.U / h
            
    #         f1x = -1/(2*height.visc) * (px_hs[i]*hx + pxx_hs[i]*h) + height.U/(h**2) * hx
    #         # v_Sj = 0 # numerical integration dy at xi
            
    #         f1x_b = -1/(3*height.visc)* h *pxx_hs[i] #eqn 10b
            
    #         f1_x_err[i] = (f1x- f1x_b)
            
    #         for j in range(height.Ny):
                
    #             y = ys[j]

    #             if y <= h:
                    
    #                 us[j,i]= 1/(2*height.visc)* pxs[j,i] * (y**2) + f1* y + height.U
    #                 vs[j,i] = -1/(6*height.visc) * pxxs[j,i] * y**3 - 1/2 * f1x_b * y**2 # eqn. 10

    #                 # ux = 1/(2*height.visc) * pxxs[j,i] * (y**2) + f1x * y 
    #                 # v_Sj = v_Sj + -ux * height.dy 
    #                 # vs[j,i] = v_Sj
                    
    #     # graphics.plot_2D(f1_x_err, height.xs, 'f1_x error', ['x', 'f1_x - f1_x'])
    #     return us, vs