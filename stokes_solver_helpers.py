# -*- coding: utf-8 -*-
"""
Created on Tue May 21 13:58:57 2024

@author: sarah
"""
import numpy as np
import graphics
# -----------------------------------------------------------------------------

# Velocity <-> Stream update

# u[k] = c2 * (psi_N - psi_S) - c3 * (u_N + u_S)
# v[k] = -c2 * (psi_E - psi_W) - c3 * (v_E + v_W)

def uv_approx(tri, u, v, psi_mirr):

    
    n = tri.Nx
    m = tri.Ny

    U = tri.U
    
    c2 = 3/(4*tri.dx)
    c3 = 1/4

    for k in range(n*m):
        i = k % n
        j = k // n

        # y=yL moving boundary
        if j == m-1: 
            u[k] = U
            v[k] = 0 
            
        # other boundaries & dead zones
        elif not tri.is_interior(i,j):
            u[k] = 0
            v[k] = 0 
                
        else: #interior
            # (u,v) at 4 point stencil

            # North (i, j+1)
            if j+1 == m-1: 
                u_N = U
                psi_N = 0
            else:
                k_N = (j+1)*n + i
                u_N = u[k_N]
                psi_N = psi_mirr[k_N]
                
            # East (i+1, j)
            if i+1 == n-1 : 
                v_E = 0
                psi_E = 0 
            else:
                k_E = j*n + i + 1
                v_E = v[k_E]
                psi_E = psi_mirr[k_E] 
            
            # South (i, j-1)
            if j-1 == 0 : 
                u_S = 0
                psi_S = 0 
            else:
                k_S = (j-1)*n + i
                u_S = u[k_S]
                psi_S = psi_mirr[k_S]
                
            # West (i-1, j)  
            if i-1 == 0 :
                v_W = 0
                psi_W = 0 
            else:
                k_W = j*n + i - 1
                v_W = v[k_W]
                psi_W = psi_mirr[k_W]
   
            u[k] = c2 * (psi_N - psi_S) - c3 * (u_N + u_S)
            v[k] = -c2 * (psi_E - psi_W) - c3 * (v_E + v_W)
    
    return u, v



# -----------------------------------------------------------------------------
# boundary mirroring

def mirror_boundary(tri, psi):
    psi_mirr = psi.copy()
    
    n = tri.Nx 
    m = tri.Ny 
    slope = tri.slope
    dx = tri.dx
    i_mid = n//2
    
    for j in range(m-1):
        
        dj = j % slope #j-index within sub triangle between two true boundary points
        di = int(j//slope) 
        
        k_left = j*n + i_mid - di
        k_right = j*n + i_mid + di
        
        if dj == 0: #true boundary point
            psi_mirr[k_left] = 0 
            psi_mirr[k_right] = 0
        
        else:
    
            # psi_mirr[k_left-1] = (1 - dj*dx/slope)*psi[k_left]
            # psi_mirr[k_right+1] = (1 + dj*dx/slope)*psi[k_right] 
            psi_mirr[k_left-1] = -psi[k_left] * dx*dj/(slope-dj)
            psi_mirr[k_right+1] = -psi[k_right] * dx*dj/(slope-dj)

    return psi_mirr
        
def unmirror_boundary(tri, psi):
    n = tri.Nx
    m = tri.Ny
    slope = tri.slope
    
    i_mid = n//2

    for j in range(m-1):
        
        di = int(j//slope) 
    
        k_left = j*n + i_mid - di
        k_right = j*n + i_mid + di
        
        psi[k_left-1] = 0
        psi[k_right+1] = 0
    
    return psi

