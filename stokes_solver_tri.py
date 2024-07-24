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
    
    M, ext_nbrs = Dpsi_cscmatrixBuild(ex)
    LU = splu(M)
    
    rhs = update_rhs(ex, u, v, old_psi, ext_nbrs)
    
    t_k0 = time.time()
    print("N=%d constr. t=%.2f"%(ex.N, t_k0-t0))
     
    for k in range(iters): 
        t_ki = time.time()
    
        psi = LU.solve(rhs)

        u, v = uv_approx(ex, u, v, psi, ext_nbrs)
        
        rhs = update_rhs(ex, u, v, psi, ext_nbrs)
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
    space = ex.space
    mat = DPsi_Mat(m, n)
    
    ext_nbrs= np.zeros((m,n,8))

    for k in range(m*n):
    
        i = k % n
        j = k // n
        
        # k = j*n + i

            
        # exterior & boundry points --> identity row
        # if not ex.is_interior(i, j):
        if space[j,i] != 1: # not interior (boundary or exterior )
            mat.append(k, k, 1)
            
        # interior 
        # [... 1 -8  1 ... -8  28 -8 ... 1 -8  1 ...]
        else: 
            mat.append(k, k, 28)
            # if a nbr is exterior... adjust on rhs each iteration
            if j+1 < m and i+1 < n and space[j+1,i+1] == -1: #NE
                ext_nbrs[j,i,0]=1
            else:
                mat.append(k, (j+1)*n + i+1, 1)
                
            if j+1 < m and space[j+1,i] == -1: #N
                ext_nbrs[j,i,1]=1
            else:
                mat.append(k, (j+1)*n + i, -8)
                
            if j+1 < m and i-1 >0 and space[j+1,i-1] == -1: #NW
                ext_nbrs[j,i,2]=1
            else:
                mat.append(k, (j+1)*n + i-1, 1)
                
            if i+1 < n and space[j,i+1] == -1: #E
                ext_nbrs[j,i,3]=1
            else:
                mat.append(k, j*n + i+1, -8)

            if i-1 > 0 and space[j,i-1] == -1: #W
                ext_nbrs[j,i,4]=1
            else:
                mat.append(k, j*n + i-1, -8)
                
            if j-1 > 0 and i+1 < n and space[j-1,i+1] == -1: #SE
                ext_nbrs[j,i,5]=1
            else:
               mat.append(k, (j-1)*n + i+1, 1)
               
            if j-1 > 0 and space[j-1,i] == -1: #S
                ext_nbrs[j,i,6]=1
            else:
                mat.append(k, (j-1)*n + i, -8)
                    
            if j-1 > 0 and i-1 > 0 and space[j-1,i-1] == -1: #SW
                ext_nbrs[j,i,7]=1
            else:
                mat.append(k, (j-1)*n + i-1, 1) 

                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           

    csc_mat = csc_matrix((mat.coefs, (mat.row, mat.col)), (m*n, m*n))

    return csc_mat, ext_nbrs


# -----------------------------------------------------------------------------
# rhs = c0*A + c1*(B - C) + Dpsi_bc

# A = u_S - u_N + v_E - v_W
# B = v_C * (u_E + u_W + u_N + u_S)
# C = u_C * (v_E + v_W + v_N + v_S)

def update_rhs(ex, u, v, psi, ext_nbrs): #
    
    n = ex.Nx 
    m = ex.Ny
    space = ex.space
    rhs = np.zeros(n*m)

    c0 = 3 * ex.dx
    c1 = 0.5 * ex.dx**2 * ex.Re
    
    for k in range(n*m):
        i = k % n
        j = k // n
        
        
        #Psi = 0 for exterior and boundary
        if space[j,i] != 1:
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
            

            dpsi_bc = 0
                                
            #possible exterior nbrs
            
            if ext_nbrs[j,i,0]: #NE:
                dpsi_bc += -1 * stream_interp_diag(psi[k])
            if ext_nbrs[j,i,1]: #N:
                dpsi_bc += 8 * stream_interp(psi[k])
            if ext_nbrs[j,i,2]: #NW:
                dpsi_bc += -1 * stream_interp_diag(psi[k])
            if ext_nbrs[j,i,3]: #E:
                dpsi_bc += 8 * stream_interp(psi[k])
            if ext_nbrs[j,i,4]: #W:
                dpsi_bc += 8 * stream_interp(psi[k])
            if ext_nbrs[j,i,5]: #SE:
                dpsi_bc += -1 * stream_interp_diag(psi[k])
            if ext_nbrs[j,i,6]: #S:
                dpsi_bc += 8 * stream_interp(psi[k])
            if ext_nbrs[j,i,7]: #SW:
                dpsi_bc += -1 * stream_interp_diag(psi[k])
                
            A = u_S - u_N + v_E - v_W
            B = v_C * (u_E + u_W + u_N + u_S)
            C = u_C * (v_E + v_W + v_N + v_S)

            rhs[k] = c0 * A + c1 * (B - C) + dpsi_bc


    return rhs
 
# Velocity <-> Stream update

# u[k] = c2 * (psi_N - psi_S) - c3 * (u_N + u_S)
# v[k] = -c2 * (psi_E - psi_W) - c3 * (v_E + v_W)

def uv_approx(tri, u, v, psi, ext_nbrs):
  

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
        elif (tri.space[j,i] != 1):
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
                if ext_nbrs[j,i,1]: #N:

                    psi_N = stream_interp(psi[k])
                else:
                    psi_N = psi[k_N]
                
            # East (i+1, j)
            if i+1 == n-1 : 
                v_E = 0
                psi_E = 0 
            else:
                k_E = j*n + i + 1
                v_E = v[k_E]
                if ext_nbrs[j,i,3]: #E:

                    psi_E = stream_interp(psi[k])

                else:
                    psi_E = psi[k_E] 
            
            # South (i, j-1)
            if j-1 == 0 : 
                u_S = 0
                psi_S = 0 
            else:
                k_S = (j-1)*n + i
                u_S = u[k_S]
                if ext_nbrs[j,i,6]: #S:

                    psi_S = stream_interp(psi[k])
                else:
                    psi_S = psi[k_S]
                
            # West (i-1, j)  
            if i-1 == 0 :
                v_W = 0
                psi_W = 0 
            else:
                k_W = j*n + i - 1
                v_W = v[k_W]
                if ext_nbrs[j,i,4]: #W:

                    psi_W = stream_interp(psi[k])

                else:
                    psi_W = psi[k_W]
   
            u[k] = c2 * (psi_N - psi_S) - c3 * (u_N + u_S)
            v[k] = -c2 * (psi_E - psi_W) - c3 * (v_E + v_W)
    
    return u, v

def stream_interp(psi_k):
    return -psi_k

def stream_interp_diag(psi_k):
    return -np.sqrt(2) * psi_k 


# s(x_a) = (psi_k) * (x_a - x_h)/(x_k-x_h)   
                                                                                                                                                                             
# s(x) = (psi_k - psi_h) * (x - x_h)/(x_k-x_h) + psi_h        
        
        
        
        
        
        
        
        
        
        