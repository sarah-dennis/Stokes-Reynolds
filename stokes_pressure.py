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
            p[k] = p[(j-1)*n+i] - py[k]* dy 
        j-=1    
       
    return p
            
    
# do them all at once : {uxx, uyy, vxx, vyy}
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
            
            if i==n-1 or space[j,i+1]==-1:# backward diff
                uxx_k = (2*u_k - 5*u[j*n+i-1] + 4*u[j*n+i-2] - u[j*n+i-3])/ex.dx**3
                vxx_k = (2*v_k - 5*v[j*n+i-1] + 4*v[j*n+i-2] - v[j*n+i-3])/ex.dx**3            
            elif i == 0 or space[j,i-1]==-1:# forward diff
                
                uxx_k = (2*u_k - 5*u[j*n+i+1] + 4*u[j*n+i+2] - u[j*n+i+3])/ex.dx**3
                vxx_k = (2*v_k - 5*v[j*n+i+1] + 4*v[j*n+i+2] - v[j*n+i+3])/ex.dx**3           
            else: #center diff
                uxx_k = (u[j*n + i+1] -2*u_k + u[j*n + i-1])/ex.dx**2
                vxx_k = (v[j*n + i+1] -2*v_k + v[j*n + i-1])/ex.dx**2

            px[k] = uxx_k + vxx_k
            
            # uyy & vyy <--| N:j+1 & S:j-1
            if j == m-1: # backward diff
                uyy_k = (2*u_k - 5*u[(j-1)*n+i] + 4*u[(j-2)*n+i] - u[(j-3)*n+i])/ex.dy**3
                vyy_k = (2*v_k - 5*v[(j-1)*n+i] + 4*v[(j-2)*n+i] - v[(j-3)*n+i])/ex.dy**3
                
            elif j == 0 or space[j-1,i]==-1: # forward diff
                uyy_k = (2*u_k - 5*u[(j+1)*n+i] + 4*u[(j+2)*n+i] - u[(j+3)*n+i])/ex.dy**3
                vyy_k = (2*v_k - 5*v[(j+1)*n+i] + 4*v[(j+2)*n+i] - v[(j+3)*n+i])/ex.dy**3

            else:#center diff
                uyy_k = (u[(j+1)*n + i] -2*u_k + u[(j-1)*n + i])/ex.dy**2
                vyy_k = (v[(j+1)*n + i] -2*v_k + v[(j-1)*n + i])/ex.dy**2

            py[k] = uyy_k + vyy_k
    return px, py












