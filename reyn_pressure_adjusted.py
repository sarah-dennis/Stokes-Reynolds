# -*- coding: utf-8 -*-
"""
Created on Tue Feb 25 15:16:58 2025

@author: sarah
"""
import numpy as np
import domain as dm
import graphics
import reyn_pressure_finDiff as fd

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
    M = fd.make_mat(height)
    rhs = adj_rhs(height, pxs, p2xs, p3xs)
    sigmas  = np.linalg.solve(M, rhs)

    U = height.U
    visc = height.visc
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

                ps_adj[j,i] = reyn_ps[i] + adj + sigmas[i]
                            

    return ps_adj, pxs, p2xs, p3xs, p4xs

def adj_rhs(height, pxs, p2xs, p3xs):
    vs = np.zeros(height.Nx)
    U = height.U
    visc = height.visc
    for i in range(height.Nx):
        h = height.hs[i]
        hx = height.hxs[i]
        h2x = height.h2xs[i]
        h3x = height.h3xs[i]
        px = pxs[i]
        
        p2x = p2xs[i]
        p3x = p3xs[i]
        
        v_a = (h**2)/2 *(px*h3x + 3*(p2x*h2x + p3x*hx) - p3x) 
        v_b = (1/2) * (px*h2x + 2*p2x*hx) * (1 + 3*hx)
        v_c = -U*visc*(h3x - h2x*(1-hx)/h - 2*hx/(h**2))
        vs[i] = -(height.dx**2)*(h**2)/visc * (v_a + v_b + v_c)
    vs[0] = 0
    vs[-1] = 0
    return vs
    