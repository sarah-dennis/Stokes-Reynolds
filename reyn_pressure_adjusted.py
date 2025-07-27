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
    
    #--------------------------------------------------------------------------
    pxs = dm.center_diff(reyn_ps, height.Nx, height.dx)
    p2xs = dm.center_second_diff(reyn_ps, height.Nx, height.dx)
    p3xs = dm.center_third_diff(reyn_ps, height.Nx, height.dx)
    p4xs = dm.center_fourth_diff(reyn_ps, height.Nx, height.dx)
    
    for i in height.i_peaks[1:-1]:
        pxs[i-1:i+2] = dm.avg_2x(pxs[i-2 : i+3])
        p2xs[i-1:i+2] = dm.avg_2x(p2xs[i-2 : i+3])
        p3xs[i-2:i+3] = dm.avg_3x(p3xs[i-3 : i+4])
        p4xs[i-3:i+4] = dm.avg_4x(p4xs[i-4 : i+5])
    
    # graphics.plot_2D_multi([pxs,p2xs,p3xs,p4xs], height.xs, 'Reynolds Pressure gradients', ['$p_x$','$p_{xx}$','$p_{xxx}$','$p_{xxxx}$'], ['x','p_{*}'])
   #---------------------------------------------------------------------------
    M = fd.make_mat(height)
    # rhs = adj_rhs(height, pxs, p2xs) #t & g velocity
    rhs = adj_rhs_new(height, pxs, p2xs, p3xs, p4xs) #new velocity
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
                adj = -pxx*(y**2)/2 - phi1x*y 

                ps_adj[j,i] = reyn_ps[i] + adj + sigmas[i]*visc 

#               ps_adj[j,i] = reyn_ps[i] + adj  #+ vy*visc

                            
    # graphics.plot_2D(sigmas, height.xs, 'sigmas', ['$x$','$\sigma(x)$'])
    return ps_adj, pxs, p2xs, p3xs, p4xs, sigmas

def adj_rhs(height, pxs, p2xs):
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

        v_a = (h**4)/2 *(3*p2x*h2x + px*h3x) 
        v_b = 3*(h**3)* (2*p2x*(hx**2) + px*h2x*hx)
        v_c = -U*visc*((h**2)*h3x - 6*(hx**3))
        vs[i] = -1/visc * (v_a + v_b + v_c)

    vs[0] = 0
    vs[-1] = 0
    return vs


def adj_rhs_new(height, pxs, p2xs, p3xs, p4xs):
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
        p4x = p4xs[i]

        v_a =(h**5)*p4x + 5*(h**4)*p3x*hx 
        v_b = (h*p4x+ 3*hx*p3x + h3x*px + 3*h2x*p2x)*(h**4) 
        v_c = 4*(h*p3x + 2*hx*p2x + h2x*px)*(h**3)*hx
        v_d = ((h**2)*h3x-2*(hx**3)-2*h*hx*h2x)
    
        vs[i] = 3*v_a/(20*visc) - (v_b + v_c)/(4*visc) + v_d*U/2

    vs[0] = 0
    vs[-1] = 0
    return vs
    