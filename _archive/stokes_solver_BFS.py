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

def run_spLU(ex, u, v, old_psi, iters, past_iters, error_mod, write_mod):
    t0 = time.time()
    
    M = Dpsi_cscmatrixBuild(ex)
    LU = splu(M)

    t_k0 = time.time()
    print("N=%d constr. t=%.2f"%(ex.N, t_k0-t0))
    
    for k in range(iters): 
        t_ki = time.time()
        
        u, v = uv_approx(ex, u, v, old_psi) 
        
        rhs = update_rhs(ex, u, v)
    
        psi = LU.solve(rhs)

        t_kj = time.time()
        
        if k % error_mod == 0: 

            max_err = np.max(np.abs(old_psi - psi))
            print(" k=%d max error: %.4e psi"%(k, max_err))
            print(" k=%d time: %.2f s"%(k, t_kj-t_ki))

        if k % write_mod == 0:
            rw.write_solution(ex, u, v, psi, k+1+past_iters)
            
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
        
        
def Dpsi_cscmatrixBuild(ex):
    m = ex.Ny
    n = ex.Nx

    mat = DPsi_Mat(m, n)

    for k in range(m*n):
    
        i = k % n
        j = k // n
        
        # k = j*n + i
        
        # exterior & boundry points --> identity row, set by rhs
        if ex.space[j,i] != 1:
    
            mat.append(k, k, 1)
            
        # interior 
        # [... 1 -8  1 ... -8  28 -8 ... 1 -8  1 ...]
        else: 
            mat.append(k, k, 28)
            
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

def update_rhs(ex, u, v): #
    
    n = ex.Nx 
    m = ex.Ny
    
    rhs = np.zeros(n*m)

    c0 = 3 * ex.dx
    c1 = 0.5 * ex.dx**2 * ex.Re
    
    for k in range(n*m):
        i = k % n
        j = k // n
        

            
        if i == 0: # inlet
            rhs[k] = ex.streamInlet(j)
        
        elif i == n-1: # outlet
            rhs[k] = ex.streamOutlet(j)
            
        elif j == m-1: # moving upper surface
            rhs[k] = ex.flux
                
        elif ex.space[j,i] != 1: #exterior to domain
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
 

def uv_approx(ex, u, v, psi):

    
    n = ex.Nx
    m = ex.Ny

    U = ex.U
    
    c2 = 3/(4*ex.dx)
    c3 = 1/4

    for k in range(n*m):
        i = k % n
        j = k // n

        # y=yL moving boundary
        if j == m-1: 
            u[k] = U
            v[k] = 0 
            
        elif i == 0: #inlet
            u[k] = ex.velInlet(j)
            v[k] = 0 
            
        elif i == n-1: #outlet

            u[k] = ex.velOutlet(j)
            v[k] = 0 
            
        # other boundaries & dead zones
        elif ex.space[j,i] != 1:
            u[k] = 0
            v[k] = 0 
                
        else: #interior
        # (u,v) at 4 point stencil
        
            if j+1 == m-1: 
                u_N = U
                psi_N = ex.flux    
            else:
                k_N = (j+1)*n + i
                u_N = u[k_N]
                psi_N = psi[k_N]

            if i+1 == n-1:
                v_E = 0 
                psi_E = ex.streamOutlet(j)
            elif ex.space[j,i+1] != 1:
                v_E = 0
                psi_E = 0
            else:
                k_E = j*n + i + 1
                v_E = v[k_E]
                psi_E = psi[k_E] 
            
            if i-1 == 0:
                v_W = 0
                psi_W = ex.streamInlet(j)
            elif ex.space[j,i-1] != 1:
                v_W = 0
                psi_W = 0
            else:
                k_W = j*n + i - 1
                v_W = v[k_W]
                psi_W = psi[k_W]
                      
            if ex.space[j-1,i] != 1:
                u_S = 0
                psi_S = 0
            else:
                k_S = (j-1)*n + i
                u_S = u[k_S]
                psi_S = psi[k_S]    
            
            u[k] = c2 * (psi_N - psi_S) - c3 * (u_N + u_S)
            v[k] = -c2 * (psi_E - psi_W) - c3 * (v_E + v_W)
        
    return u, v        
        
        

        
        
        
        
        
        