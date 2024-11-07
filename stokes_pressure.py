# -*- coding: utf-8 -*-
"""
Created on Wed Aug 28 11:58:25 2024

@author: sarah
"""
import numpy as np

def resistance(ex, p):
    j_in = int((ex.yf - ex.H_in/2)/ex.dy)
    j_out =int((ex.yf - ex.H_out/2)/ex.dy)
    
    dp = p[j_out*ex.Nx+ex.Nx-1] - p[j_in*ex.Nx+1]

    if ex.flux!=0:
        r= dp/ex.flux
        return dp, r
    else:
        return dp, np.inf


def pressure(ex, u, v):
    px, py = fishfun(ex, u, v)
    n = ex.Nx
    m = ex.Ny
    dy = ex.dy
    dx = ex.dx
    shape = n*m
    p = np.zeros(shape)

    # set the ambient pressure at outlet 
    p[(m-2)*n + n-1] = ex.p_ambient    
    # contour the first interior row (backwards) using px
    i=n-2
    while i >= 0:
        k =   (m-2)*n + i
        k_E = (m-2)*n + i+1
        p[k] = p[k_E]-px[k]*dx
        i-=1
     
    # set the top boundary using dy from interior row
    # p[(m-1)*n-1] = ex.p_ambient
    i = n-1
    while i >=0:
        k   = (m-1)*n + i
        k_S = (m-2)*n + i
        p[k] = p[k_S] + py[k_S]*dy
        i-=1

    for i in range(n):
        j=m-3
        while j >= 0:
            k = j*n + i
            if ex.space[j,i]==1:
                k_N = (j+1)*n + i
                p[k] = p[k_N] - py[k_N]*dy
            j-=1    
    
    for i in range(n):
        j=m-3
        while j >= 0:
            k = j*n + i
            if ex.space[j,i]==0:
                if ex.space[j+1,i]==1:
                    k_N = (j+1)*n + i
                    p[k] = p[k_N] - py[k_N]*dy
                    
                elif i < n-1 and ex.space[j,i+1] == 1:
                    k_E = j*n +i+1
                    p[k] = p[k_E] - px[k_E]*dx
                    
                elif i > 0 and ex.space[j,i-1]==1:
                    k_W = j*n+i-1
                    p[k] = p[k_W] + px[k_W]*dx
                    
                elif i < n-1 and ex.space[j+1,i+1]==1:
                    k_NE = (j+1)*n+i+1
                    p[k] = p[k_NE] - px[k_NE]*dx -py[k_NE]*dy
                    
                elif i > 0 and ex.space[j+1,i-1]==1:
                    k_NW = (j+1)*n+i-1
                    p[k] = p[k_NW] + px[k_NW]*dx -py[k_NW]*dy
                else:        
                    assert False, "Grid misalignment - check x_peaks for dx=1/N"
                
            j-=1    
    
    # note: if stream-velocity has not converged, pressure will not be continuous
    return p
            
    
# {uxx, uyy, vxx, vyy}
def fishfun(ex, u, v):
# u[k] = u[jn + i] = u(xi, yj)
    n = ex.Nx
    m = ex.Ny
    shape = n*m
    space = ex.space #[0: boundry, 1: interior, -1: exterior]
    
    px = np.zeros(shape)
    py = np.zeros(shape)
    for k in range(shape):
        i = k % n
        j = k // n
        
        if space[j,i] != 1: # exterior & boundary
            continue
        else:
            
            u_k = u[k]
            v_k = v[k]
            k_W=j*n + i-1
            k_E=j*n + i+1
            # uxx & vxx <--| E:i+1 & W:i-1 
            if space[j,i+1]==-1:
                u_E = ex.interp_E(i,j, u[k_W])
                v_E = ex.interp_E(i,j, v[k_W])
            else:
                u_E = u[k_E]
                v_E = v[k_E]
                
            if space[j,i-1] ==-1:
                u_W = ex.interp_W(i,j,  u[k_E])
                v_W = ex.interp_W(i,j,  v[k_E])
            else:
                u_W = u[k_W]
                v_W = v[k_W]
                
                
            uxx_k = (u_E -2*u_k + u_W)/ex.dx**2
            vxx_k = (v_E -2*v_k + v_W)/ex.dx**2
            
            ux_k = (u_E - u_W)/(2*ex.dx)
            vx_k = (v_E - v_W)/(2*ex.dx)


            # uyy & vyy <--| N:j+1 & S:j-1
            k_N=(j+1)*n + i
            k_S=(j-1)*n + i
            
            u_N = u[k_N]                    
            v_N = v[k_N]
            
            if space[j-1,i] == -1:
                u_S = ex.interp_S(i,j, u[k_N])
                v_S = ex.interp_S(i,j, v[k_N])
            else: 
                u_S = u[k_S]
                v_S = v[k_S]

            uyy_k = (u_N -2*u_k + u_S)/ex.dx**2
            vyy_k = (v_N -2*v_k + v_S)/ex.dx**2
            
            uy_k= (u_N - u_S)/(2*ex.dx)
            vy_k= (v_N - v_S)/(2*ex.dx)

            px[k] = ex.visc*(uxx_k + uyy_k) - ex.dens*(u_k*ux_k + v_k*uy_k)
            py[k] = ex.visc*(vxx_k + vyy_k) - ex.dens*(u_k*vx_k + v_k*vy_k)
            
    return px, py

