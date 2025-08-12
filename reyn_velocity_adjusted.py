# -*- coding: utf-8 -*-
"""
Created on Thu May  8 14:11:23 2025

@author: sarah
"""
import numpy as np
import graphics
import domain as dm


def make_adj_velocity(height, BC, adj_pressure):

    # derivatives of the reynolds pressure
    pxs = adj_pressure.reyn_pxs
    p2xs = adj_pressure.reyn_p2xs
    p3xs = adj_pressure.reyn_p3xs
    p4xs = adj_pressure.reyn_p4xs
    
    # sigmas = adj_pressure.sigmas
    sigma_xs = adj_pressure.sigma_xs
    sigma_2xs = adj_pressure.sigma_2xs

# 

    us = make_us(height, BC.U, pxs, p2xs, p3xs, sigma_xs)
    vs = make_vs(height,  BC.U, pxs, p2xs, p3xs, p4xs, sigma_xs, sigma_2xs)
      
    return us, vs


def make_us(height, U, pxs, p2xs, p3xs, sigmaxs):
    us = np.zeros((height.Ny, height.Nx))

    visc = 1#height.visc
    for i in range(height.Nx):
        h = height.hs[i]
        hx = height.hxs[i]
        h2x = height.h2xs[i]
        px = pxs[i]
        p2x = p2xs[i]
        p3x = p3xs[i]
        sigmax = sigmaxs[i]
        
        h2 = h**2
        h3 = h**3
        hx2 = hx**2
        for j in range(height.Ny):
            y = height.ys[j]
            y2 = y**2
            y3 = y**3
            y4 = y**4
            
            if y > h:
                continue
            else:
              # reynolds u
                uA = 1/(2*visc)  * px                   * (y2 - h*y) + U*(1 - y/h)
              # adjusted u
                uB = -1/(24*visc)* p3x                  * (y4 - h3*y)
                uC = 1/(12*visc) * (h*p3x + 2*p2x*hx + px*h2x)  * (y3 - h2*y)
                uD = -U/6        * (-2*hx2/h3 + h2x/h2) * (y3 - h2*y)
                uE = (1/2)       * sigmax               * (y2 - h*y)
                
                us[j,i] = (uA+uB+uC+uD+uE)

    return us


def make_vs(height, U, pxs, p2xs, p3xs, p4xs, sigmaxs, sigma2xs):
    vs = np.zeros((height.Ny, height.Nx))
    visc = 1# height.visc
    for i in range(height.Nx):
        h = height.hs[i]
        hx = height.hxs[i]
        h2x = height.h2xs[i]
        h3x = height.h3xs[i]
        px = pxs[i]
        p2x = p2xs[i]
        p3x = p3xs[i]
        p4x =  p4xs[i]
        sigmax = sigmaxs[i]
        sigma2x = sigma2xs[i]
        
        h2 = h**2
        h3 = h**3
        h4 = h**4
        hx3 = hx**3
        for j in range(height.Ny):
            y = height.ys[j]
            y2 = y**2
            y3 = y**3
            y4 = y**4
            y5 = y**5
            if y > h:
                continue
            else:
                
                # v(x,y) = integral -du/dx  dy
                vA = -1/(2*visc)* (p2x*(y3/3-h*y2/2) - px*hx*y2/2) - U*hx*y2/(2*h2) 
                vB = 1/(24*visc)*(p4x*(y5/5-h3*y2/2) - 3*h2*y2*p3x*hx/2 )
                vCa = (h*p4x + 3*p3x*hx + 3*p2x*h2x + px*h3x)*(y4/4 - h2*y2/2)
                vCb = -(h*p3x + 2*p2x*hx + px*h2x)*h*hx*y2
                vC = -1/(12*visc)*(vCa+vCb)
                vDa = (y4-2*h2*y2)/(4*h2)*h3x + (3*y4-2*h2*y2)/(2*h4)*hx3 + (4*h2*y2-3*y4)/(2*h3)*hx*h2x

                vD = U/6 * vDa
                vE = -1/2 * (sigma2x*(y3/3 -h*y2/2) - sigmax*hx*y2/2)
                vs[j,i] = vA+vB+vC+vD+vE

    return vs

    
# # ------------------------------------------------------------------------------

# as in Takeuchi-Gu            
def make_adj_velocity_TG(height, BC, adj_pressure):
    ps = adj_pressure.ps_2D #ps = reyn_ps + adj_ps
    
    pxs, p2xs = make_px_pxx(height, ps)
    
    px_hs, p2x_hs, pxy_hs = make_pxh_pxxh_pxyh(height, pxs, p2xs)
        
    us, vs = make_us_vs(height, BC, pxs, px_hs, p2xs, p2x_hs, pxy_hs)
    return us, vs

def make_px_pxx(height, ps):  
    ys = height.ys
    hs = height.hs
    dx = height.dx
    dy = height.dy
    pxs = np.zeros((height.Ny, height.Nx))
    pxxs = np.zeros((height.Ny, height.Nx))
    
    # make px and pxx
    
    for i in range(height.Nx):
        h = hs[i]
        for j in range(height.Ny):
            y = ys[j]

            if y <= h:
    
                if i < height.Nx-1 and y <= hs[i+1] and i > 0 and y <= hs[i-1]: # all interior nbrs
                    px = dm.center_first(dx, ps[j,i-1 : i+2])
                    pxx = dm.center_second(dx, ps[j,i-1 : i+2]) 
                
                elif i < height.Nx-1 and y <= hs[i+1] : # West out of bounds, fwd diff 
                    
                    if i < height.Nx-2 and y <= hs[i+2]: # 2 nbrs East
                        px = dm.right_first(dx, ps[j,i : i+3]) 
                        
                        if i < height.Nx-3 and y <= hs[i+3]: # 3 nbrs East
                            pxx = dm.right_second(dx, ps[j,i : i+4])

                        else: # only 2 nbrs East
                            pxx = dm.right_second_O1(dx, ps[j,i : i+3]) 

                    else: # only 1 nbrs East
                        px = dm.right_first_O1(dx, ps[j,i: i+2]) 
                        pxx = 0
                
                elif i > 0 and y <= hs[i-1]: # East out of bounds, bkwd diff 
                    
                    if i > 1 and y <= hs[i-2]: # 2 nbrs West
                        px = dm.left_first(dx, ps[j,i-2 : i+1]) 
                        
                        if i > 2 and y <= hs[i-3]: # 3 nbrs West
                            pxx = dm.left_second(dx, ps[j,i-3 : i+1])
                            
                        else: # only 2 nbrs West
                            pxx = dm.left_second_O1(dx, ps[j,i-2 : i+1])
                            
                    else: # only 1 nbr West
                        px = dm.left_first_O1(dx, ps[j,i-1 : i+1]) 
                        pxx = 0
                
                else: # both East and West out of bounds
                    pxy =  (pxs[j-1,i] - pxs[j-3,i])/(2*dy)
                    px = pxs[j-1,i] + pxy*dy
                    
                    pxxy =  (pxxs[j-1,i] - pxxs[j-3,i])/(2*dy)
                    pxx = pxxs[j-1,i] + pxxy*dy
            
            else: # exterior
                continue
            
            
            pxs[j,i] = px
            pxxs[j,i] = pxx
               
    # discontinuity averaging 
    for i in height.i_peaks[1:-1]:
        for j in range(height.Ny):
            h = height.ys[j]
            pxs[j,i-2: i+3] = dm.avg_3x(pxs[j,i-3 : i+4])
            pxxs[j,i-3 : i+4] = dm.avg_4x(pxxs[j,i-4 : i+5])
   
    return pxs, pxxs     

def make_pxh_pxxh_pxyh(height,pxs,pxxs):
    dy = height.dy
    
    px_hs = np.zeros(height.Nx)
    pxx_hs = np.zeros(height.Nx)
    pxy_hs = np.zeros(height.Nx)    
    
    for i in range(height.Nx):
        h = height.hs[i]
        for j in range(height.Ny):
            y = height.ys[j]
            
            if y == h :
                px_hs[i] = pxs[j,i]
                pxx_hs[i] = pxxs[j,i]
                pxy_hs[i] = dm.left_first(dy, pxs[j-2 : j+1, i])
            elif (y < h and y+dy > h):
                px_hs[i] = pxs[j,i] + pxxs[j,i]*(h-y)
                pxxy_i = dm.left_first(dy, pxxs[j-2 : j+1, i])
                pxx_hs[i] = pxxs[j,i] + pxxy_i*(h-y)
                pxy_hs[i] = dm.left_first(dy, pxs[j-2 : j+1, i]) + pxxy_i*(h-y)

    return px_hs, pxx_hs, pxy_hs

def make_us_vs(height, BC, pxs, px_hs, p2xs, p2x_hs, pxy_hs):
    
    us = np.zeros((height.Ny, height.Nx))
    vs = np.zeros((height.Ny, height.Nx))

    U = BC.U
    visc = 1# height.visc
    for i in range(height.Nx):
        h = height.hs[i]
        hx = height.hxs[i]
        px_h = px_hs[i]
        pxx_h = p2x_hs[i]
        pxy_h = pxy_hs[i]
        
        phi1  = -1/(2*visc)*h*px_h - U/h
        
        phi1x = -1/(2*visc)*(h*(pxx_h + pxy_h*hx) + hx*px_h) + U/(h**2) *hx

        # when using sigmas = 0 in p_adj, v_h corrects the velocity boundary
        v_h = 1/(12*visc)*(h**3) * (pxx_h+ 3*pxy_h*hx  + 3*px_h*hx/h) -U/2 * hx
        
        for j in range(height.Ny):
            y = height.ys[j]
            px = pxs[j,i]
            pxx = p2xs[j,i]
            
            if y < h:
                us[j,i]= 1/(2*visc)*px* y**2  + phi1* y + U
                
                vs[j,i] = -1/(6*visc)*pxx* y**3 - 1/2*phi1x* y**2 -v_h*y/h            

            else:
                continue

    return us, vs


    