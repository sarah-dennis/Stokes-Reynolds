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
    
        for i in range(3, self.height.Nx-2):
            qs[i]= np.sum(vx[:,i])*self.height.dx

        # graphics.plot_2D(qs[3:-3], self.height.xs[3:-3], 'flux', ['x', 'q'])
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
                    
                    if y <= h_E and y <= h_W:  # interior

                        ux = domain.center_first(dx, u[j,i-1 : i+2]) #(p_E - p_W)/(2*height.dx)
                        
                    elif y <= h_E and y > h_W: # West out of bounds, fwd diff (right sided)

                        if y <= h_EE:
                            ux = domain.right_first(dx, u[j, i : i+3]) #(-3*p +4*p_E -p_EE)/(2*height.dx)

                        else:
                            ux = domain.right_first_O1(dx, u[j, i: i+2]) #(p_E - p)/height.dx

                    else: # East out of bounds, bkwd diff (left sided)

                        if y <= h_WW:
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
        # graphics.plot_contour_mesh(inc, height.xs, height.ys, 'incompressibility', ['I', 'x', 'y'], -1, 1)
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
        dx = height.dx
        hs = height.hs
        ys = height.ys
        us = np.zeros((height.Ny, height.Nx))
        vs = np.zeros((height.Ny, height.Nx))
        pxs = np.zeros((height.Ny, height.Nx))
        pxxs = np.zeros((height.Ny, height.Nx))
        px_hs = np.zeros(height.Nx)
        pxx_hs = np.zeros(height.Nx)
    
    
        for j in range(height.Ny):
           
            y = ys[j]
            pj = ps[j]
            
            for i in range(3,height.Nx-3):
                 
               
                h = hs[i]

                if y <= h:
                    if y <= hs[i+1] and y <= hs[i-1]: #interior nbrs
                        px = domain.center_first(dx, pj[i-1 : i+2])
                        pxx = domain.center_second(dx, pj[i-1 : i+2]) 
                    
                    elif y <= hs[i+1] and y > hs[i-1]: # West out of bounds, fwd diff (right sided)
                        if y <= hs[i+2]:
                            px = domain.right_first(dx, pj[i : i+3]) 
                            
                            if y <= hs[i+3]:
                                pxx = domain.right_second(dx, pj[i : i+4])
                                
                            else:
                                pxx = domain.right_second_O1(dx, pj[i : i+3]) 
    
                        else:
                            px = domain.right_first_O1(dx, pj[i: i+2]) 
                            pxx = 0
                    
                    elif y > hs[i+1] and y <= hs[i-1]: # East out of bounds, bkwd diff (left sided)
                        if y <= hs[i-2]:
                            px = domain.left_first(dx, pj[i-2 : i+1]) 
                            
                            if y <= hs[i-3]:
                                pxx = domain.left_second(dx, pj[i-3 : i+1]) 
                            else:
                                pxx = domain.left_second_O1(dx, pj[i-2 : i+1])
                                
                        else:
                            px = domain.left_first_O1(dx, pj[i-1 : i+1]) 
                            pxx = 0
                           
                    else: # both East and West out of bounds
                        px = 0
                        pxx = 0
                
                else: #exterior
                    px = 0
                    pxx = 0
                
                pxs[j,i] = px
                pxxs[j,i] = pxx
                
                if y == h or (y < h and ys[j+1] >= h):
                    px_hs[i] = px
                    pxx_hs[i] = pxx

        # print(px_hs, pxx_hs)
        # graphics.plot_contour_mesh(pxs, height.xs, height.ys, 'pxs', ['$p_{x}$', '$x$', '$y$'], vmin=-3, vmax=3)
        # graphics.plot_contour_mesh(pxxs, height.xs, height.ys, 'pxxs', ['$p_{xx}$', '$x$', '$y$'], vmin=-3, vmax=3)

        for i in range(height.Nx):
            h = hs[i]
            # hx = height.hxs[i]
            
            
            f1 =  -1/(2*height.visc) * h * px_hs[i] - height.U / h
            
            # f1x = -1/(2*height.visc) * (px_hs[i]*hx + pxx_hs[i]*h) + height.U/(h**2) * hx
            # v_Sj = 0 # numerical integration dy at xi
            
            f1x = -1/(3*height.visc)* h *pxx_hs[i] #eqn 10b
            
            for j in range(height.Ny):
                
                y = ys[j]

                if y <= h:
                    
                    us[j,i]= 1/(2*height.visc)* pxs[j,i] * (y**2) + f1* y + height.U
                    vs[j,i] = -1/(6*height.visc) * pxxs[j,i] * y**3 - 1/2 * f1x * y**2 # eqn. 10

                    # ux = 1/(2*height.visc) * pxxs[j,i] * (y**2) + f1x * y 
                    # v_Sj = v_Sj + -ux * height.dy 
                    # vs[j,i] = v_Sj
                
                    
                    
        # graphics.plot_contour_mesh(us, height.xs, height.ys, '$u(x,y)$', ['$u$', '$x$', '$y$'], vmin=-1, vmax=1)
        # graphics.plot_contour_mesh(vs, height.xs, height.ys, '$v(x,y)$', ['$v$', '$x$', '$y$'], vmin=-1, vmax=1)

        uv_mag = np.sqrt(us**2 + vs**2)
        graphics.plot_contour_mesh(uv_mag, height.xs, height.ys, '$|(u,v)|_2$', ['$|(u,v)|_2$', '$x$', '$y$'], vmin=0, vmax=3, log_cmap=False)

        # vy=np.flip(vy, 0)
        # vx=np.flip(vx, 0)

        return us, vs