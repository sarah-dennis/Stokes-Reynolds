# -*- coding: utf-8 -*-
"""
Created on Wed Aug 28 11:58:25 2024

@author: sarah
"""
import numpy as np

def get_dp(ex, p):
    
    p_2D = p.reshape((ex.Ny,ex.Nx))


    dp_dim = (sum(p_2D[:,0]) - sum(p_2D[:,-1]))*ex.dy

    return dp_dim


def pressure(ex, u, v):
    px, py = px_py(ex, u, v)    

    n = ex.Nx
    m = ex.Ny

    shape = n*m
    p = np.zeros(shape)
    
    
    L_x = (ex.xf-ex.x0)
    dy = ex.dy/L_x
    dx = ex.dx/L_x
    if ex.flux == 0:
        U_star = ex.U
    else:   
        U_star = ex.flux/L_x
    P_star = U_star/L_x
    
    
    p[2*n-1] = ex.p_ambient  # i=n-1, j=1
    #k = j*n + i
    #k = n + (n-1)
    
    i=n-2 #j=1
    while i >= 0:
        k = n + i
        k_E = n + i+1
        p[k] = p[k_E] - px[k]*dx
        i-=1
     
    i=n-2 #j=0
    while i >= 0:
        k = i
        k_N = n + i
        p[k] = p[k_N] - py[k_N]*dy
        i-=1
            
    for i in range(n):
        j = 2
        while j <=m-1 and ex.space[j,i] != -1:
            k = j*n + i
            
            k_S = (j-1)*n + i
            p[k] = p[k_S] + py[k_S]*dy
          
            j+=1  
    p=p*P_star
    return p
            
    
# {uxx, uyy, vxx, vyy,}
def px_py(ex, u, v):
# u[k] = u[jn + i] = u(xi, yj)
    n = ex.Nx
    m = ex.Ny
    shape = n*m
    space = ex.space #[0: boundry, 1: interior, -1: exterior]
    
    L_x = (ex.xf-ex.x0)
    if ex.flux == 0:
        U_star = ex.U
    else:   
        U_star = ex.flux/L_x
    u = u/U_star
    v = v/U_star
    dx = ex.dx/L_x
    
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
            
            if space[j,i+1]==-1:
                scale_E = ex.scale_E(i,j)
                u_E = ex.interp(scale_E, u[k_W])
                v_E = ex.interp(scale_E, v[k_W])
            else:
                u_E = u[k_E]
                v_E = v[k_E]
                
            if space[j,i-1] ==-1:
                scale_W = ex.scale_W(i,j)
                u_W = ex.interp(scale_W, u[k_E])
                v_W = ex.interp(scale_W, v[k_E])
            else:
                u_W = u[k_W]
                v_W = v[k_W]
                 
            uxx_k = (u_E -2*u_k + u_W)/dx**2
            vxx_k = (v_E -2*v_k + v_W)/dx**2
            
            ux_k = (u_E - u_W)/(2*dx)
            vx_k = (v_E - v_W)/(2*dx)

            # uyy & vyy <--| N:j+1 & S:j-1
            k_N=(j+1)*n + i
            k_S=(j-1)*n + i
            
            u_S = u[k_S]
            v_S = v[k_S]
            
            if space[j+1,i] == -1:
                scale_N = ex.scale_N(i,j)
                u_N = ex.interp(scale_N, u[k_S])
                v_N = ex.interp(scale_N, v[k_S])
                
            else: 
                u_N = u[k_N]
                v_N = v[k_N]
            
            uyy_k = (u_N -2*u_k + u_S)/dx**2
            vyy_k = (v_N -2*v_k + v_S)/dx**2
            
            uy_k = (u_N - u_S)/(2*dx)
            vy_k = (v_N - v_S)/(2*dx)

            
            px[k] = (uxx_k + uyy_k) - ex.Re* (u_k*ux_k + v_k*uy_k)
            py[k] = (vxx_k + vyy_k) - ex.Re* (u_k*vx_k + v_k*vy_k)

    return px, py

