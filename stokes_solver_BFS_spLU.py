# -*- coding: utf-8 -*-
"""
Created on Fri May 24 08:49:45 2024

@author: sarah
"""
import numpy as np
import time
import stokes_readwrite as rw
from scipy.sparse import csc_matrix
from scipy.sparse.linalg import splu

def run_spLU(tri, u, v, old_psi, iters, past_iters, error_mod, write_mod):
    t0 = time.time()
    
    M = Dpsi_cscmatrixBuild(tri)
    LU = splu(M)

    t_k0 = time.time()
    print("N=%d constr. t=%.2f"%(tri.N, t_k0-t0))
    
    for k in range(iters): 
        t_ki = time.time()
        
        u, v = uv_approx(tri, u, v, old_psi) 
        
        rhs = update_rhs(tri, u, v)
    
        psi = LU.solve(rhs)

        t_kj = time.time()
        
        if k % error_mod == 0: 

            max_err = np.max(np.abs(old_psi - psi))
            print(" k=%d max error: %.4e psi"%(k, max_err))
            print(" k=%d time: %.2f s"%(k, t_kj-t_ki))

        if k % write_mod == 0:
            rw.write_solution(tri, u, v, psi, k+1+past_iters)
            
        old_psi = psi

    return u, v, psi

class DPsi_Mat():
    def __init__(self, Ny, Nx):
        self.length = Ny * Nx * 9
        self.coefs = np.zeros(self.length)
        self.row = np.zeros(self.length)
        self.col = np.zeros(self.length)
        self.s = 0
        
    def append(self, k, nbr_k, coef):
        self.row[self.s] = k
        self.col[self.s] = nbr_k
        self.coefs[self.s] = coef
        self.s += 1
        
        
def Dpsi_cscmatrixBuild(tri):
    m = tri.Ny
    n = tri.Nx

    mat = DPsi_Mat(m, n)

    for k in range(m*n):
    
        i = k % n
        j = k // n
        
        # k = j*n + i
        
        # exterior & boundry points --> identity row, set by rhs
        if not tri.is_interior(i, j):
    
            mat.append(k, k, 1)
            
        # interior 
        # [... 1 -8  1 ... -8  28 -8 ... 1 -8  1 ...]
        else: 
            mat.append(k, k, 28)
            
            # if a nbr is on the boundary... nothing changes, the boundary is solved with identity rows
            # with rect-lin shape nbrs don't fall in the exterior
            
            mat.append(k, j*n + i-1, -8)
            mat.append(k, j*n + i+1, -8)
            mat.append(k, (j-1)*n + i, -8)
            mat.append(k, (j+1)*n + i, -8)
                
            mat.append(k, (j-1)*n + i-1, 1)  
            mat.append(k, (j-1)*n + i+1, 1)
            mat.append(k, (j+1)*n + i-1, 1)
            mat.append(k, (j+1)*n + i+1, 1)
            
    csc_mat = csc_matrix((mat.coefs, (mat.row, mat.col)), (m*n, m*n))

    return csc_mat


# -----------------------------------------------------------------------------
# rhs = c0*A + c1*(B - C) 
# A = u_S - u_N + v_E - v_W
# B = v_C * (u_E + u_W + u_N + u_S)
# C = u_C * (v_E + v_W + v_N + v_S)

def update_rhs(tri, u, v): #
    
    n = tri.Nx 
    m = tri.Ny
    
    rhs = np.zeros(n*m)

    c0 = 3 * tri.dx
    c1 = 0.5 * tri.dx**2 * tri.Re
    
    for k in range(n*m):
        i = k % n
        j = k // n
        
        if j == m-1: # moving upper surface
            rhs[k]= tri.flux
            
        elif i == 0: # inlet
            rhs[k] = tri.streamInlet(j)
        
        elif i == n-1: # outlet
            rhs[k] = tri.streamOutlet(j)

        elif not tri.is_interior(i,j):
            rhs[k] = 0
        
        else: # interior
            # (u,v) at 9 point stencil

            u_C = u[k]
            v_C = v[k]

            # North (i, j+1) 
            k_N = (j+1)*n + i
            u_N = u[k_N]
            v_N = v[k_N]
                
            # East (i+1, j)                
            k_E = j*n + i + 1
            u_E = u[k_E]
            v_E = v[k_E]
                
            # South (i, j-1)
            k_S = (j-1)*n + i
            u_S = u[k_S]
            v_S = v[k_S]
                
            # West (i-1, j)          
            k_W = j*n + i - 1
            u_W = u[k_W]
            v_W = v[k_W]
        
    
            A = u_S - u_N + v_E - v_W
            B = v_C * (u_E + u_W + u_N + u_S)
            C = u_C * (v_E + v_W + v_N + v_S)

            rhs[k] = c0 * A + c1 * (B - C)



    return rhs
 

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
            
        elif i == 0: #inlet
            u[k] = tri.velInlet(j)
            v[k] = 0 
            
        elif i == n-1: #outlet

            u[k] = tri.velOutlet(j)
            v[k] = 0 
            
        # other boundaries & dead zones
        elif not tri.is_interior(i,j):
            u[k] = 0
            v[k] = 0 
                
        else: #interior
            # (u,v) at 4 point stencil
            if j+1 == m-1:
                u_N = U
                psi_N = tri.flux
            else:
                k_N = (j+1)*n + i
                u_N = u[k_N]
                psi_N = psi[k_N]
          
            if tri.is_lowerbndry(i,j):
                u_S = 0
                psi_S = 0
            else:
                k_S = (j-1)*n + i
                u_S = u[k_S]
                psi_S = psi[k_S]    
            
            if i+1 == n-1:
                v_E = 0 
                psi_E = tri.streamOutlet(j)
            else:
                k_E = j*n + i + 1
                v_E = v[k_E]
                psi_E = psi[k_E] 
            
            if i-1 == 0:
                v_W = 0
                psi_W = tri.streamInlet(j)
            else:
                k_W = j*n + i - 1
                v_W = v[k_W]
                psi_W = psi[k_W]

            u[k] = c2 * (psi_N - psi_S) - c3 * (u_N + u_S)
            v[k] = -c2 * (psi_E - psi_W) - c3 * (v_E + v_W)
    
    return u, v        
        
        

        
        
        
        
        
        