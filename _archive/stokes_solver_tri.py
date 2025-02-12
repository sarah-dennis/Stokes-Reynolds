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

def run_spLU(ex, u, v, old_psi, iters, past_iters, error_mod, write_mod, err_tol):
    t0 = time.time() 
    M = Dpsi_cscmatrixBuild(ex)
    LU = splu(M)
    
    t_k0 = time.time()
    print("N=%d constr. t=%.2f"%(ex.N, t_k0-t0))
     
    for k in range(iters): 
        t_ki = time.time()
        
        u, v = uv_approx(ex, u, v, old_psi)
        
        rhs = update_rhs(ex, u, v, old_psi)
    
        psi = LU.solve(rhs)

        t_kj = time.time()
        
        if k % error_mod == 0: 

            max_err = np.max(np.abs(old_psi - psi))
            print(" k=%d max error: %.4e psi"%(k, max_err))
            print(" k=%d time: %.2f s"%(k, t_kj-t_ki))
            if max_err < err_tol:
                break
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

    for k in range(m*n):
    
        i = k % n
        j = k // n
        
        # k = j*n + i 

            
        # exterior -> identity row = 0
        # boundary -> identiy row = b.c.
        if space[j,i] != 1:
            mat.append(k, k, 1)
            
        # interior 
        # [... 1 -8  1 ... -8  28 -8 ... 1 -8  1 ...]
        else: 
            mat.append(k, k, 28)
            # if a nbr is exterior |-> adjust on rhs each iteration
            # if j+1 < m and i+1 < n and space[j+1,i+1] == -1:
            #     ext_nbrs[j,i,0]=1
            # else:
            mat.append(k, (j+1)*n + i+1, 1) #NE
                
            # if j+1 < m and space[j+1,i] == -1:
            #     ext_nbrs[j,i,1]=1
            # else:
            mat.append(k, (j+1)*n + i, -8) #N
                
            # if j+1 < m and i-1 > 0 and space[j+1,i-1] == -1:
            #     ext_nbrs[j,i,2]=1
            # else:
            mat.append(k, (j+1)*n + i-1, 1) #NW
                
            # if i+1 < n and space[j,i+1] == -1:
            #     ext_nbrs[j,i,3]=1
            # else:
            mat.append(k, j*n + i+1, -8) #E

            # if i-1 > 0 and space[j,i-1] == -1: 
            #     ext_nbrs[j,i,4]=1
            # else:
            mat.append(k, j*n + i-1, -8)#W
                
            # if j-1 > 0 and i+1 < n and space[j-1,i+1] == -1: 
            #     ext_nbrs[j,i,5]=1
            # else:
            mat.append(k, (j-1)*n + i+1, 1)#SE
               
            # if j-1 > 0 and space[j-1,i] == -1:
            #     ext_nbrs[j,i,6]=1
            # else:
            mat.append(k, (j-1)*n + i, -8) #S
                    
            # if j-1 > 0 and i-1 > 0 and space[j-1,i-1] == -1: 
            #     ext_nbrs[j,i,7]=1
            # else:
            mat.append(k, (j-1)*n + i-1, 1) #SW

                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           
    csc_mat = csc_matrix((mat.coefs, (mat.row, mat.col)), (m*n, m*n))

    return csc_mat


# -----------------------------------------------------------------------------
# rhs = c0*A + c1*(B - C) + Dpsi_bc

# A = u_S - u_N + v_E - v_W
# B = v_C * (u_E + u_W + u_N + u_S)
# C = u_C * (v_E + v_W + v_N + v_S)

def update_rhs(ex, u, v, psi): #
    
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
        if j == m-1:
            rhs[k] = ex.flux
            
        elif space[j,i] != 1:
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
            

            
            psi_k = psi[k]          
            
            #possible exterior nbrs 
            dpsi_bc = 0

            if ex.space[j+1,i+1] == -1: #NE:
                dpsi_bc += -1 * ex.stream_interp_EW(j+1, i, psi[(j+1)*n+i])

            if ex.space[j+1,i-1] == -1: #NW:
                dpsi_bc += -1 * ex.stream_interp_EW(j+1, i, psi[(j+1)*n+i])
                                                     
            if ex.space[j,i+1] == -1: #E:
                dpsi_bc += 8 * ex.stream_interp_EW(j, i, psi_k)
                
            if ex.space[j,i-1] == -1: #W:
                dpsi_bc += 8 * ex.stream_interp_EW(j, i, psi_k)
                
            if ex.space[j-1,i] == -1: #S:
                dpsi_bc += 8 * ex.stream_interp_S(j,i, psi_k)
                
            if ex.space[j-1,i+1] == -1: #SE:
                
                if i+1 < ex.apex:
                    slope = ex.slope_A
                else:
                    slope = ex.slope_B
                
                if abs(slope) < 1:
                    
                    dpsi_bc += -1 * ex.stream_interp_S(j, i+1, psi[j*n+i+1])
                else:
                    dpsi_bc += -1 * ex.stream_interp_EW(j-1, i, psi[(j-1)*n+i]) 

            if ex.space[j-1,i-1]==-1: #SW:
                                
                if i-1 < ex.apex:
                    slope = ex.slope_A
                else:
                    slope = ex.slope_B
                
                if abs(slope) < 1:
                    
                    dpsi_bc += -1 * ex.stream_interp_S(j, i-1, psi[j*n+i-1])
                else:
                    dpsi_bc += -1 * ex.stream_interp_EW(j-1, i, psi[(j-1)*n+i]) 
                    
                    
            A = u_S - u_N + v_E - v_W
            B = v_C * (u_E + u_W + u_N + u_S)
            C = u_C * (v_E + v_W + v_N + v_S)

            rhs[k] = c0 * A + c1 * (B - C) + dpsi_bc


    return rhs
 
# Velocity <-> Stream update

# u[k] = c2 * (psi_N - psi_S) - c3 * (u_N + u_S)
# v[k] = -c2 * (psi_E - psi_W) - c3 * (v_E + v_W)

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
            
        # other boundaries & exterior
        elif (ex.space[j,i] != 1):
            u[k] = 0
            v[k] = 0 
                
        else: #interior
            # (u,v, psi) at 4 point stencil
            
            # North (i, j+1)
            if j+1 == m-1: 
                u_N = U
                psi_N = ex.flux
            else:
                k_N = (j+1)*n + i
                u_N = u[k_N]
                psi_N = psi[k_N]
                
            # East (i+1, j)
            if ex.space[j,i+1] == -1: 
                v_E = 0
                psi_E = ex.stream_interp_EW(j,i,psi[k])
            elif ex.space[j,i+1] == 0:
                v_E = 0
                psi_E = 0
            else:
                k_E = j*n + i + 1
                v_E = v[k_E]
                psi_E = psi[k_E] 

            # West (i-1, j)  
            if ex.space[j,i-1] == -1:
                v_W = 0
                psi_W = ex.stream_interp_EW(j,i,psi[k])
            elif ex.space[j,i-1] == 0:
                v_W = 0
                psi_W = 0 
            else:
                k_W = j*n + i - 1
                v_W = v[k_W]
                psi_W = psi[k_W]
   
            # South (i, j-1)
            if ex.space[j-1,i] == -1:
                u_S = 0 
                psi_S = ex.stream_interp_S(j,i,psi[k])
            elif ex.space[j-1,i] == 0:
                u_S = 0
                psi_S = 0 
            else:
                k_S = (j-1)*n + i
                u_S = u[k_S]
                psi_S = psi[k_S]
                

            u[k] = c2 * (psi_N - psi_S) - c3 * (u_N + u_S)
            v[k] = -c2 * (psi_E - psi_W) - c3 * (v_E + v_W)
    
    return u, v


        
        
        
        
        
        
        
        
        