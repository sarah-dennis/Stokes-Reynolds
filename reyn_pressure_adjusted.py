# -*- coding: utf-8 -*-
"""
Created on Tue Feb 25 15:16:58 2025

@author: sarah
"""
import numpy as np
import domain as dm
import graphics


def make_adj_ps(height, reyn_ps):
    ps_adj = np.zeros((height.Ny, height.Nx))

    errs = np.zeros(height.Nx)
    
    hs = height.hs
    hxs = height.hxs
    dx = height.dx
    dy = height.dy
    #--------------------------------------------------------------------------
    pxs = np.zeros(height.Nx)
    pxxs = np.zeros(height.Nx)

    
    i=0
    pxs[i] =  dm.right_first(dx, reyn_ps[i : i+3])
    pxxs[i] = dm.right_second(dx, reyn_ps[i : i+4])

 
    i=height.Nx-1
    pxs[i] = dm.left_first(dx, reyn_ps[i-2 : i+1])
    pxxs[i] = dm.left_second(dx, reyn_ps[i-3 : i+1])

    
    for i in range(1,height.Nx-1): 
        pxs[i] = dm.center_first(dx, reyn_ps[i-1 : i+2])
        pxxs[i] = dm.center_second(dx, reyn_ps[i-1 : i+2])

    for i in height.i_peaks[1:-1]:
        pxs[i-1:i+2] = dm.avg_2x(pxs[i-2 : i+3])
        pxxs[i-1:i+2] = dm.avg_2x(pxxs[i-2 : i+3])
   #---------------------------------------------------------------------------
    visc = height.visc
    U = height.U
    for i in range(height.Nx): 
        h = hs[i]
        hx = hxs[i]
        px = pxs[i]
        pxx = pxxs[i]
        
        phi1x = -(pxx*h + px*hx)/2 + U*visc/(h**2)*hx
                
        vy = (h/(2*visc)*px - U/h)*hx
        
        for j in range(height.Ny):
            y = height.ys[j]
            if y > h:
                ps_adj[j,i] = None
                continue
            else:  
                adj = -pxx*(y**2)/2 - phi1x*y - vy*visc

                ps_adj[j,i] = reyn_ps[i] + adj
                            
                if y + height.dy > h:
                    
                    errs[i] = abs(adj)
    print('max err p(x,h)=p_re(x): ', np.max(errs))
    return ps_adj 