# -*- coding: utf-8 -*-
"""
Created on Thu May  8 14:11:23 2025

@author: sarah
"""
import numpy as np
import graphics
import domain
        
# def make_adj_velocity(height, ps):
    
#     # v(x,y) = (int_0^y[p] dy - int_0^h [p] dy * y)/visc + v(x,0) + vy(x,0)y
#     vs= make_vs(height, ps)

#     # u(x,y) = int_x0^x[du/dx] dx
#     us = make_us(height, vs, ps)
   
#     return us, vs


# def make_vs(height, ps):
#     vs = np.zeros((height.Ny, height.Nx))
#     intPdys, phis = make_intPdys(height, ps)
    
#     for i in range(height.Nx):
#         h = height.hs[i]
        
#         for j in range(height.Ny):
#             y = height.ys[j]
#             if y <= h:
#                 if y == 0:
#                     vs[j,i] = 0 
#                 else:
#                     vs[j,i] =  1/height.visc * (intPdys[j,i] -ps[0,i]*y - phis[i]*y)
#             else:
#                 vs[j,i] = 0
#                 continue

#     return vs

# # for v(x,y)

# def make_intPdys(height, ps):
#     intPdys = np.zeros((height.Ny, height.Nx)) # int_y0^y[p] dy
#     phis = np.zeros( height.Nx) # int_y0^h[p] dy

#     for i in range(height.Nx):
#         h = height.hs[i]
#         intPdy_i =0
        
#         for j in range(1,height.Ny):
#             y = height.ys[j]
    
#             if y <= h:
#                 intPdy_i += (ps[j-1,i] + ps[j,i]) *height.dy/2          
#                 intPdys[j,i] = intPdy_i
            
#             # integration constant: int_0^h p(x,y) dy
#             if y == h:
#                 phis[i] = intPdy_i/h - ps[0,i]
#                 continue
#             elif y < h and y + height.dy >= h : # boundary approximation 
#                 phis[i]= intPdy_i/h - ps[0,i] + ps[j,i]*(h-y)/h
#                 continue
   
    
#     return intPdys, phis


# def make_us(height, vs, ps):
#     vys = make_dvdys(height, vs)

#     h0 = height.hs[0]
#     px0 = domain.right_first(height.dx, ps[0,0:3])
#     q = (height.U*h0)/2 - (px0*(h0**3))/(12*height.visc)

#     us = np.zeros((height.Ny, height.Nx))
#     pwl_reg = 0 # x=h(y)
    
#     for i in range(height.Nx):
#         h = height.hs[i]
#         x = height.xs[i]
        
#         if i == height.i_peaks[pwl_reg] and i < height.Nx-1:# skip i=0~x0
#             pwl_reg +=1
            
#         for j in range(height.Ny):
#             y = height.ys[j]

#             if j == 0: # y=y0
#                 us[j,i] = height.U
        
#             elif y == h:
#                 us[j,i] = 0 
            
#             #vertical boundary
#             elif i == height.i_peaks[pwl_reg-1]  and y >= height.hs[i] and y <= height.hs[i+1]:
#                 us[j,i] = 0
            
#             elif y < h:
            
#                 if i == 0: #inlet
#                     us[j,i] = height.U*(h-y)*(h-3*y)/h**2 + 6*q*y*(h-y)/h**3
                    
#                 elif y <= height.hs[i-1]: 
                    
#                     if i > 1 and y <= height.hs[i-2]: # contour using i-1 and i-2
#                         ux = -vys[j,i]
#                         us[j,i] = (4*us[j,i-1] -us[j,i-2] +2*ux*height.dx)/3 
                        
#                     else: # contour using i-1 only
#                         ux =  -vys[j,i]
#                         us[j,i] = us[j,i-1] +ux * height.dx
                
                                    
#                 elif i-1 == height.i_peaks[pwl_reg-1]: # i-1 vertical boundary, contour using i-1 only
#                     ux = - vys[j,i]
#                     us[j,i] =  us[j,i-1] +ux * height.dx 
                       
#                 else: # (i-1,j) out of bounds, interpolate
#                     ux = - vys[j,i]
#                     if i < height.Nx-1 and h - height.hs[i+1] != 0:
#                         x_bdry_W = height.xs[i+1] + height.dx * (y - height.hs[i+1])/(h - height.hs[i+1])
#                     else:
#                         x_bdry_W = x    
#                     us[j,i] = ux * (x-x_bdry_W)
                    
#             else:
#                 us[j,i] = 0

            
#     return us

# # for u(x,y)...          
# def make_dvdys(height, vs):
#     vys = np.zeros((height.Ny,height.Nx))

#     for i in range(1,height.Nx):
#         h = height.hs[i] 
#         vys[0,i]= domain.right_first(height.dy,vs[0:3,i]) 

#         for j in range(1,height.Ny):
#             y = height.ys[j]    
#             if y <= h:
#                 if j < height.Ny-1 and height.ys[j+1] <= h:
#                     vys[j,i] = domain.center_first(height.dy,vs[j-1:j+2,i]) 
                    
                    
#                 else:
                    
#                     vys[j,i] = domain.left_first(height.dy,vs[j-2:j+1,i]) 
                    
#             else:
#                 vys[j,i] = 0

#     return vys

#------------------------------------------------------------------------------
# as in Takeuchi-Gu            
def make_adj_velocity_old(height, ps):
    
    pxs, pxxs = make_px_pxx(height, ps)

    px_hs, pxx_hs, pxy_hs = make_pxh_pxxh_pxyh(height, pxs, pxxs)
    
    us, vs = make_us_vs(height, pxs, pxxs, px_hs, pxx_hs, pxy_hs)
    return us, vs

def make_px_pxx(height, ps):  
    ys = height.ys
    hs = height.hs
    dx = height.dx
    
    pxs = np.zeros((height.Ny, height.Nx))
    pxxs = np.zeros((height.Ny, height.Nx))
    
    # make px and pxx
    for j in range(height.Ny):
        y = ys[j]
        for i in range(height.Nx):
            h = hs[i]

            if y <= h:
    
                if i < height.Nx-1 and y <= hs[i+1] and i > 0 and y <= hs[i-1]: # all interior nbrs
                    px = domain.center_first(dx, ps[j,i-1 : i+2])
                    pxx = domain.center_second(dx, ps[j,i-1 : i+2]) 
                
                elif i < height.Nx-1 and y <= hs[i+1] : # West out of bounds, fwd diff 
                
                    if i < height.Nx-2 and y <= hs[i+2]: # 2 nbrs East
                        px = domain.right_first(dx, ps[j,i : i+3]) 
                        
                        if i < height.Nx-3 and y <= hs[i+3]: # 3 nbrs East
                            pxx = domain.right_second(dx, ps[j,i : i+4])

                        else: # only 2 nbrs East
                            pxx = domain.right_second_O1(dx, ps[j,i : i+3]) 

                    else: # only 1 nbrs East
                        px = domain.right_first_O1(dx, ps[j,i: i+2]) 
                        pxx = 0
                
                elif i > 0 and y <= hs[i-1]: # East out of bounds, bkwd diff 
                
                    if i > 1 and y <= hs[i-2]: # 2 nbrs West
                        px = domain.left_first(dx, ps[j,i-2 : i+1]) 
                        
                        if i > 2 and y <= hs[i-3]: # 3 nbrs West
                            pxx = domain.left_second(dx, ps[j,i-3 : i+1])
                            
                        else: # only 2 nbrs West
                            pxx = domain.left_second_O1(dx, ps[j,i-2 : i+1])
                            
                    else: # only 1 nbr West
                        px = domain.left_first_O1(dx, ps[j,i-1 : i+1]) 
                        pxx = 0
                
                else: # both East and West out of bounds
                    px = pxs[j-1,i] + pxxs[j-1,i]*dx
                    pxx = pxxs[j-1,i]
            
            else: # exterior
                px = 0
                pxx = 0
            
            pxs[j,i] = px
            pxxs[j,i] = pxx
        
    # discontinuity averaging
    for i in height.i_peaks[1:-1]:
        for j in range(height.Ny):
            pxs[j,i-2 : i+3] = domain.avg_3x(pxs[j,i-3 : i+4])
            pxxs[j,i-5 : i+6] = domain.avg_6x(pxxs[j,i-6 : i+7])
    return pxs, pxxs     

def make_pxh_pxxh_pxyh(height,pxs,pxxs):
    dy = height.dy
    
    px_hs = np.zeros(height.Nx)
    pxx_hs = np.zeros(height.Nx)
    pxy_hs = np.zeros(height.Nx)    
    #integration constants px(x,h) and pxx(x,h) 
    for i in range(height.Nx):
        h = height.hs[i]
        for j in range(height.Ny):
            y = height.ys[j]
            
            
            if y == h :
                px_hs[i] = pxs[j,i]
                pxx_hs[i] = pxxs[j,i]
                pxy_hs[i] = domain.left_first(dy, pxs[j-2 : j+1, i])
            elif (y < h and y+dy > h):
                px_hs[i] = pxs[j,i] + pxxs[j,i]*(h-y)
                pxxy_i = domain.left_first(dy, pxxs[j-2 : j+1, i])
                pxx_hs[i] = pxxs[j,i] + pxxy_i*(h-y)
                pxy_hs[i] = domain.left_first(dy, pxs[j-2 : j+1, i]) + pxxy_i*(h-y)
                
    return px_hs, pxx_hs, pxy_hs

def make_us_vs(height, pxs, pxxs, px_hs, pxx_hs, pxy_hs):
    
    
    us = np.zeros((height.Ny, height.Nx))
    vs = np.zeros((height.Ny, height.Nx))

    # u and v
    # max_re=0
    max_err=0
    U = height.U
    visc = height.visc
    for i in range(height.Nx):
        h = height.hs[i]
        hx = height.hxs[i]
        px_h = px_hs[i]
        pxx_h = pxx_hs[i]
        pxy_h = pxy_hs[i]
        
        phi1  = -1/(2*visc)*h*px_h - U/h
        
        phi1x = -1/(2*visc)*(h*(pxx_h + pxy_h*hx) + hx*px_h) + U/(h**2) *hx

        
        # f1x = -1/(6*height.visc)*(h**2)*pxx_hs[i] -h/2 * phi1x
        re = 1/(12*visc)*(h**3) * (pxx_h+ 3*pxy_h*hx  + 3*px_h*hx/h) -height.U/2 * hx


        
        for j in range(height.Ny):
            y = height.ys[j]
            px = pxs[j,i]
            pxx = pxxs[j,i]
            
            if y < h:
                us[j,i]= 1/(2*visc)*px* y**2  + phi1* y + U
                
                vs[j,i] = -1/(6*visc)*pxx* y**3 - 1/2*phi1x* y**2 #-re*y/h            
                if y + height.dy :
                    
                    if abs(vs[j,i])> max_err:
                        max_err = abs(vs[j,i])
            else:
                continue
            
    print('max err v(x,h)=0', max_err)
    return us, vs

    # for i in range(height.Nx):
    #     h = hs[i]
    #     hx = height.hxs[i]
    #     f1  = -1/(2*height.visc)*h*pxs[0,i] - height.U / h
    #     f1x = -1/(2*height.visc)*(h*pxxs[0,i] +hx*pxs[0,i]) + height.U/(h**2) *hx
  
    #     for j in range(height.Ny):
    #         y = ys[j]
    #         if y <= h:
    #             us[j,i]= 1/(2*height.visc)*pxs[0,i]* y**2  + f1* y + height.U
    #             vs[j,i] = -1/(6*height.visc)*pxxs[0,i]* y**3 - 1/2*f1x* y**2 
    #         else:
    #             continue

    