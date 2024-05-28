# -*- coding: utf-8 -*-
"""
Created on Fri May 24 08:49:45 2024

@author: sarah
"""
import numpy as np
import time
import stokes_readwrite as rw
from stokes_solver_helpers import uv_approx, mirror_boundary
from scipy.sparse import csc_matrix
from scipy.sparse.linalg import splu

#TODO: changed order so mirror then uv approx
def run_spLU(tri, u, v, old_psi, iters, past_iters, error_mod, write_mod):
    t0 = time.time()
    M = Dpsi_cscmatrixBuild(tri)
    LU = splu(M)
    psi_mirr = mirror_boundary(tri, old_psi)
    rhs = update_rhs(tri, u, v, psi_mirr)
    t_k0 = time.time()
    print("N=%d constr. t=%.2f"%(tri.N, t_k0-t0))
     
    for k in range(iters): 
        t_ki = time.time()
    
        psi = LU.solve(rhs)

        psi_mirr = mirror_boundary(tri, psi)

        u, v = uv_approx(tri, u, v, psi_mirr)
        
        rhs = update_rhs(tri, u, v, psi_mirr)
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

    # slope = tri.slope
    # nc = n//2
    mat = DPsi_Mat(m, n)
    # coefs = np.zeros(m*n*9)
    # row = np.zeros(m*n*9)
    # col = np.zeros(m*n*9)

    
    # s=0
    
    for k in range(m*n):
    
        i = k % n
        j = k // n
        
        # k = j*n + i
        
        # Identity row for exterior points
        if not tri.is_interior(i, j):
            mat.append(k, k, 1)
            # row[s] = k
            # col[s] = k
            # coefs[s] = 1
            # s+=1
        
        # [... 1 -8  1 ... -8  28 -8 ... 1 -8  1 ...]
        else: #i,j is interior --> check 8 nbrs 
            mat.append(k, k, 28)
            # row[s] = k
            # col[s] = k
            # coefs[s] = 28
            # s+=1
            
            if tri.is_interior(i-1, j):
                mat.append(k, j*n + i-1, -8)
                # row[s] = k
                # col[s] = j*n + i - 1
                # coefs[s] = -8
                # s+=1 
            
            if tri.is_interior(i+1, j):
                mat.append(k, j*n + i+1, -8)
                # row[s] = k
                # col[s] = j*n + i + 1
                # coefs[s] = -8
                # s+=1 
            
            if tri.is_interior(i, j-1):
                mat.append(k, (j-1)*n + i, -8)
                # row[s] = k
                # col[s] = (j-1)*n + i
                # coefs[s] = -8
                # s+=1
            
            
            if tri.is_interior(i, j+1):
                mat.append(k, (j+1)*n + i, -8)
                # row[s] = k
                # col[s] = (j+1)*n + i
                # coefs[s] = -8
                # s+=1
                
            if tri.is_interior(i-1, j-1):
                mat.append(k, (j-1)*n + i-1, 1)  
                # row[s] = k
                # col[s] = (j-1)*n + i-1
                # coefs[s] = 1
                # s+=1

            if tri.is_interior(i+1, j-1):
                mat.append(k, (j-1)*n + i+1, 1)
                # row[s] = k
                # col[s] = (j-1)*n + i+1
                # coefs[s] = 1
                # s+=1

            if tri.is_interior(i-1, j+1):
                mat.append(k, (j+1)*n + i-1, 1)
                # row[s] = k
                # col[s] = (j+1)*n + i-1
                # coefs[s] = 1
                # s+=1

            if tri.is_interior(i+1, j+1):
                mat.append(k, (j+1)*n + i+1, 1)
                # row[s] = k
                # col[s] = (j+1)*n + i+1
                # coefs[s] = 1
                # s+=1
            
            # if nbr is exterior... adjust on rhs (each iteration, relies on updated psi)

    csc_mat = csc_matrix((mat.coefs, (mat.row, mat.col)), (m*n, m*n))

    return csc_mat


# -----------------------------------------------------------------------------
# rhs = c0*A + c1*(B - C) + Dpsi_bc

# A = u_S - u_N + v_E - v_W
# B = v_C * (u_E + u_W + u_N + u_S)
# C = u_C * (v_E + v_W + v_N + v_S)

def update_rhs(tri, u, v, psi_mirr): #
    
    n = tri.Nx 
    # nc = tri.apex # i-index of triangle apex
    m = tri.Ny
    
    rhs = np.zeros(n*m)


    # slope = tri.slope
    U = tri.U
    
    c0 = 3 * tri.dx
    c1 = 0.5 * tri.dx**2 * tri.Re
    
    for k in range(n*m):
        i = k % n
        j = k // n

        # boundary & exterior zones
        # if i == 0 or j == 0 or i == n-1 or j == m-1:
        #     rhs[k] = 0
               
        # elif (i < nc and j < (m-1) - slope*i) or (i > nc and j < slope * (i - nc)):
        #     rhs[k] = 0
        
        if not tri.is_interior(i,j):
            rhs[k] = 0
        
        # interior
        else:

            # (u,v) at 9 point stencil
            #k = j*n + i
            u_C = u[k]
            v_C = v[k]
            dpsi_bc = 0
            
            # North (i, j+1) -- never exterior
            k_N = (j+1)*n + i
            if j+1 == m-1: 
                u_N = U
                v_N = 0
            else:
                u_N = u[k_N]
                v_N = v[k_N]
                
                
            # East (i+1, j)                
            k_E = j*n + i + 1
            if not tri.is_interior(i+1,j): 
                u_E = 0
                v_E = 0
                dpsi_bc += 8 * psi_mirr[k_E]
            else:
                u_E = u[k_E]
                v_E = v[k_E]
                
            # South (i, j-1)
            k_S = (j-1)*n + i
            if not tri.is_interior(i, j-1): 
                u_S = 0
                v_S = 0
                dpsi_bc += 8 * psi_mirr[k_S]
            else:
                u_S = u[k_S]
                v_S = v[k_S]
                
            # West (i-1, j)          
            k_W = j*n + i - 1
            if not tri.is_interior(i-1, j): 
                u_W = 0
                v_W = 0
                dpsi_bc += 8 * psi_mirr[k_W]
            else:

                u_W = u[k_W]
                v_W = v[k_W]
                
            #NorthEast
            k_NE = (j+1)*n + i+1
            if not tri.is_interior(i+1, j+1): 
                dpsi_bc += -1 * psi_mirr[k_NE]

            #NorthWest
            k_NW = (j+1)*n + i-1
            if not tri.is_interior(i-1, j+1): 
                dpsi_bc += -1 * psi_mirr[k_NW]

            #SouthEast
            k_SE = (j-1)*n + i+1
            if not tri.is_interior(i+1, j-1): 
                dpsi_bc += -1 * psi_mirr[k_SE]

            #SouthWest
            k_SW = (j-1)*n + i-1
            if not tri.is_interior(i-1, j-1): 
                dpsi_bc += -1 * psi_mirr[k_SW]

            A = u_S - u_N + v_E - v_W
            B = v_C * (u_E + u_W + u_N + u_S)
            C = u_C * (v_E + v_W + v_N + v_S)

            rhs[k] = c0 * A + c1 * (B - C) + dpsi_bc

    return rhs
 

        
        
        
        
        
        
        
        
        
        
        
        