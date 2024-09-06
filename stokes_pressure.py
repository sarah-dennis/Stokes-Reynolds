# -*- coding: utf-8 -*-
"""
Created on Wed Aug 28 11:58:25 2024

@author: sarah
"""
import numpy as np

def vorticity(ex, u_2D, v_2D):
    uy_2D = np.gradient(u_2D, ex.dy, axis=0)
    vx_2D = np.gradient(v_2D, ex.dx, axis=1)
    w_2D = np.zeros((ex.Ny,ex.Nx))
    for j in range(ex.Ny):
        for i in range(ex.Nx):   
            w_2D[j,i] = vx_2D[j,i] - uy_2D[j,i]
    return w_2D

def pressure(ex, u, v):
    px, py = fishfun(ex, u, v)
    n = ex.Nx
    m = ex.Ny
    dy = ex.dy
    dx = ex.dx
    shape = n*m
    p = np.zeros(shape)
    
    k0 = (m-1)*n
    p[k0] = ex.p_ambient
    
    for i in range(1,n):
        k=(m-1)*n + i
        p[k] = p[k-1]+ px[k]* dx
        
    j= m-2
    while j >=0:
        for i in range(n):
            k=j*n + i
            p[k] = p[(j+1)*n+i] - py[k]*dy 
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
        
        if space[j,i] == -1: # exterior
            continue
        else:
            
            u_k = u[k]
            v_k = v[k]
    
            # uxx & vxx <--| E:i+1 & W:i-1 
            if i == 0 or space[j,i-1]==-1:
                u_E = u[j*n + i+1]
                u_2E = u[j*n + i+2]
                u_3E = u[j*n + i+3]
                
                v_E = v[j*n + i+1]
                v_2E = v[j*n + i+2]
                v_3E = v[j*n + i+3]
                
                uxx_k = (2*u_k - 5*u_E + 4*u_2E - u_3E)/ex.dx**3
                vxx_k = (2*v_k - 5*v_E + 4*v_2E - v_3E)/ex.dx**3
                
            elif i == n-1 or space[j,i+1] ==-1:
                u_W = u[j*n + i-1]
                u_2W = u[j*n + i-2]
                u_3W = u[j*n + i-3]
                
                v_W = v[j*n + i-1]
                v_2W = v[j*n + i-2]
                v_3W = v[j*n + i-3]
                
                uxx_k = (2*u_k - 5*u_W + 4*u_2W - u_3W)/ex.dx**3
                vxx_k = (2*v_k - 5*v_W + 4*v_2W - v_3W)/ex.dx**3
                
            else:
            

                u_W = u[j*n + i-1]
                v_W = v[j*n + i-1]
                
                u_E = u[j*n + i+1]
                v_E = v[j*n + i+1]
                
                uxx_k = (u_E -2*u_k + u_W)/ex.dx**2
                vxx_k = (v_E -2*v_k + v_W)/ex.dx**2


            
            # uyy & vyy <--| N:j+1 & S:j-1
            

            if j==m-1:
                u_S = u[(j-1)*n + i]
                u_2S = u[(j-2)*n + i]
                u_3S = u[(j-3)*n + i]
                
                v_S = v[(j-1)*n + i]
                v_2S = v[(j-2)*n + i]
                v_3S = v[(j-3)*n + i]
                
                uyy_k = (2*u_k - 5*u_S + 4*u_2S - u_3S)/ex.dy**3
                vyy_k = (2*v_k - 5*v_S + 4*v_2S - v_3S)/ex.dy**3
                
            elif j==0 or space[j-1,i] == -1:
                u_N = u[(j+1)*n + i]
                u_2N = u[(j+2)*n + i]
                u_3N = u[(j+3)*n + i]
                
                v_N = v[(j+1)*n + i]
                v_2N = v[(j+2)*n + i]
                v_3N = v[(j+3)*n + i]
                
                uyy_k = (2*u_k - 5*u_N + 4*u_2N - u_3N)/ex.dy**3
                vyy_k = (2*v_k - 5*v_N + 4*v_2N - v_3N)/ex.dy**3
                    
            else: 

                u_S = u[(j-1)*n + i]
                v_S = v[(j-1)*n + i]
                
                u_N = u[(j+1)*n + i]
                v_N = v[(j+1)*n + i]
                
                uyy_k = (u_N -2*u_k + u_S)/ex.dy**2
                vyy_k = (v_N -2*v_k + v_S)/ex.dy**2
                
            px[k] = uxx_k + uyy_k
            py[k] = vxx_k + vyy_k
    return px, py












