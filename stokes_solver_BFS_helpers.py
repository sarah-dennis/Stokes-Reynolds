# -*- coding: utf-8 -*-
"""
Created on Tue May 21 13:58:57 2024

@author: sarah
"""

# -----------------------------------------------------------------------------

# Velocity <-> Stream update

# u[k] = c2 * (psi_N - psi_S) - c3 * (u_N + u_S)
# v[k] = -c2 * (psi_E - psi_W) - c3 * (v_E + v_W)

def uv_approx(tri, u, v, psi):

    
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
            
        elif i == 0:
            y = j*n + i
            u[k] = tri.velInlet(y)
            v[k] = 0 
        elif i == n-1:
            y = j*n + i
            u[k] = tri.velOutlet(y)
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
                psi_N = 0# tri.flux
            else:
                k_N = (j+1)*n + i
                u_N = u[k_N]
                psi_N = psi[k_N]
                
            # East (i+1, j)
            if i+1 == n-1 : 
                v_E = 0
                y = j*n + i+1
                # psi_E = 0
                tri.streamOutlet(y)
            else:
                k_E = j*n + i + 1
                v_E = v[k_E]
                psi_E = psi[k_E] 
            
            # West (i-1, j)  
            if i-1 == 0 :
                v_W = 0
                y = j*n + i-1
                # psi_W = 0
                tri.streamInlet(y)
                
            else:
                k_W = j*n + i - 1
                v_W = v[k_W]
                psi_W = psi[k_W]
   
            # South (i, j-1)
            if j-1 == 0 : 
                u_S = 0
                psi_S = 0
                
            else:
                k_S = (j-1)*n + i
                u_S = u[k_S]
                psi_S = psi[k_S]
                

            u[k] = c2 * (psi_N - psi_S) - c3 * (u_N + u_S)
            v[k] = -c2 * (psi_E - psi_W) - c3 * (v_E + v_W)
    
    return u, v

