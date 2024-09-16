# -*- coding: utf-8 -*-
"""
Created on Wed Aug 28 11:58:25 2024

@author: sarah
"""
import numpy as np

#TODO: bcs??
def vorticity(ex, u_2D, v_2D):
    uy_2D = np.gradient(u_2D, ex.dy, axis=0)
    vx_2D = np.gradient(v_2D, ex.dx, axis=1)
    w_2D = np.zeros((ex.Ny,ex.Nx))
    for j in range(ex.Ny):
        for i in range(ex.Nx):   
            w_2D[j,i] = vx_2D[j,i] - uy_2D[j,i]
    return w_2D

def resistance(ex, p):
    j_in = int((ex.yf - ex.H_in/2)/ex.dy)
    j_out =int((ex.yf - ex.H_out/2)/ex.dy)
    
    dp = p[j_out*ex.Nx+ex.Nx-1] - p[j_in+ex.Nx]
    
    if ex.flux!=0:
        R = dp/ex.flux
        return R
    else:
        return dp

def pressure(ex, u, v):
    px, py = fishfun(ex, u, v)
    n = ex.Nx
    m = ex.Ny
    dy = ex.dy
    dx = ex.dx
    shape = n*m
    p = np.zeros(shape)

    # set the ambient pressure at first interior pint (i=0
    p[(m-2)*n] = ex.p_ambient    
    
    # contour the first interior row using px
    for i in range(1,n):
        k =   (m-2)*n + i
        k_W = (m-2)*n + i-1
        p[k] = p[k_W]+ px[k]*dx
     
    # set the top boundary using dy from interior row
    p[(m-1)*n] = ex.p_ambient
    for i in range(1,n):
        k   = (m-1)*n + i
        k_S = (m-2)*n + i
        p[k] = p[k_S] + py[k_S]*dy
    

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
                    print(i,j)
            j-=1    
    
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
    
            # uxx & vxx <--| E:i+1 & W:i-1 
            if space[j,i+1]==-1:
                u_E = ex.interp_E_W(i,j, i+1,j, u_k)
                v_E = ex.interp_E_W(i,j, i+1,j, v_k)
            else:
                u_E = u[j*n + i+1]
                v_E = v[j*n + i+1]
                
            if space[j,i-1] ==-1:
                u_W = ex.interp_E_W(i,j, i-1,j, u_k)
                v_W = ex.interp_E_W(i,j, i-1,j, v_k)
            else:
                u_W = u[j*n + i-1]
                v_W = v[j*n + i-1]
                
                
            uxx_k = (u_E -2*u_k + u_W)/ex.dx**2
            vxx_k = (v_E -2*v_k + v_W)/ex.dx**2


            
            # uyy & vyy <--| N:j+1 & S:j-1

            u_N = u[(j+1)*n + i]                    
            v_N = v[(j+1)*n + i]
            
            if space[j-1,i] == -1:
                u_S = ex.interp_S(i,j, i,j-1,u_k)
                v_S = ex.interp_S(i,j, i,j-1,v_k)
            else: 
                u_S = u[(j-1)*n + i]
                v_S = v[(j-1)*n + i]

            uyy_k = (u_N -2*u_k + u_S)/ex.dy**2
            vyy_k = (v_N -2*v_k + v_S)/ex.dy**2
                
            px[k] = uxx_k + uyy_k
            py[k] = vxx_k + vyy_k
            
    return px, py












