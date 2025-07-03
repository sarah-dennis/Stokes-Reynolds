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

    hs = height.hs
    hxs = height.hxs
    dx = height.dx
    dy = height.dy
    #--------------------------------------------------------------------------
    pxs = dm.center_diff(reyn_ps, height.Nx, height.dx)
    p2xs = dm.center_second_diff(reyn_ps, height.Nx, height.dx)
    p3xs = dm.center_third_diff(reyn_ps, height.Nx, height.dx)
    p4xs = dm.center_fourth_diff(reyn_ps, height.Nx, height.dx)
  
    for i in height.i_peaks[1:-1]:
        pxs[i-1:i+2] = dm.avg_2x(pxs[i-2 : i+3])
        p2xs[i-1:i+2] = dm.avg_2x(p2xs[i-2 : i+3])
        p3xs[i-2:i+3] = dm.avg_4x(p3xs[i-3 : i+4])
        p4xs[i-2:i+3] = dm.avg_4x(p4xs[i-3 : i+4])
   #---------------------------------------------------------------------------
    visc = height.visc
    U = height.U
    for i in range(height.Nx): 
        h = hs[i]
        hx = hxs[i]
        px = pxs[i]
        pxx = p2xs[i]
        
        phi1x = -(pxx*h + px*hx)/2 + U*visc/(h**2)*hx
            
        # vy = -(px*h/(2*visc) - U/h)*hx       

        for j in range(height.Ny):
            y = height.ys[j]
            if y > h:
                ps_adj[j,i] = None
                continue
            else:  
                adj = -pxx*(y**2)/2 - phi1x*y #+ vy*visc

                ps_adj[j,i] = reyn_ps[i] + adj
                            

    return ps_adj, pxs, p2xs, p3xs, p4xs
