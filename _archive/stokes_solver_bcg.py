# -*- coding: utf-8 -*-
"""
Created on Fri May 24 08:48:38 2024

@author: sarah
"""
import numpy as np
import time

import stokes_readwrite as rw
from stokes_solver_helpers import uv_approx, mirror_boundary, unmirror_boundary
from scipy.sparse.linalg import LinearOperator
from scipy.sparse.linalg import bicgstab

def run_bicgstab(tri, u, v, past_psi, iters, past_iters, rtol, error_mod, write_mod):
    
    M = Dpsi_linOp(tri)
    
    for i in range(iters): 
        
        t0 = time.time()
        rhs = update_rhs(tri, u, v)

        psi, exit_flag = bicgstab(M, rhs, tol=rtol)
        
        u, v = uv_approx(tri, u, v, psi)
        
        psi = mirror_boundary(tri, psi)

        tf = time.time()
        if i % error_mod == 0: 
            past_psi = unmirror_boundary(tri, past_psi)
            psi = unmirror_boundary(tri, psi)
            
            err_i = np.max(np.abs(psi - past_psi))
            
            print("k=%d of %d"%(i+past_iters+1, iters+past_iters))
            print("  time: %.3f s"%(tf-t0))
            print("  error: %.5e psi"%err_i)
            
            if err_i < rtol:
                print('warning: lower bicgstab_rtol')
            
        if i % write_mod == 0:
            psi = unmirror_boundary(tri, psi)
            rw.write_solution(tri.filename, tri.Nx*tri.Ny, u, v, psi, i+1+past_iters)

        past_psi = psi

    return u, v, psi

# -----------------------------------------------------------------------------
# LinOp dPsi : 9 point discr. 
# -----------------------------------------------------------------------------

class Dpsi_linOp(LinearOperator):

    def __init__(self, tri):
        self.tri = tri
        nm = tri.Nx * tri.Ny
        self.shape = ((nm, nm))        
        self.dtype = np.dtype('f8')
        self.mv = np.zeros(nm) 
        
    def _matvec(self, psi): # v:= psi[j*n + i]  
        n = self.tri.Nx
        nc = n//2
        m = self.tri.Ny
        
        slope = self.tri.slope

        for k in range(n*m):
            i = k % n
            j = k // n
            
            # cavity boundary & dead zones
            if i == 0 or i == n-1 or j == 0 or j == m-1: 
                self.mv[k] = 0
            
            elif  (i < nc and j < (m-1) - slope*i) or (i > nc and j < slope * (i - nc)):
                self.mv[k] = 0
       
            #interior   
            else: 
                # psi[k] at 9 point stencil         
                
                psi_C = psi[k]
        
                #North (i, j+1)
                if j+1 == m-1: 
                    psi_N = 0
                else:
                    psi_N = psi[(j+1)*n + i]
                    
    
                #South (i, j-1)
                if j-1 == 0 : 
                    psi_S = 0
                else:
                    psi_S = psi[(j-1)*n + i]
                    
                #East (i+1, j)
                if i+1 == n-1:
                    psi_E = 0
                else:
                    psi_E = psi[j*n + i+1]
                    
                #West (i-1,j)
                if i-1 == 0 :
                    psi_W = 0
                else:
                    psi_W = psi[j*n + i-1]
                    
                #NorthEast (i+1, j+1)
                if i+1 == n-1 or j+1 == m-1 : 
                    psi_NE = 0
                else:
                    psi_NE  = psi[(j+1)*n + i+1]
                
                #NorthWest (i-1, j+1)
                if i-1 == 0 or j+1 == m-1: 
                    psi_NW = 0
                else:
                    psi_NW = psi[(j+1)*n + i-1]
                
                #SouthEast (i+1, j-1)
                if i+1 == n-1 or j-1 == 0: 
                    psi_SE = 0
                else:
                    psi_SE = psi[(j-1)*n + i+1]
                
                #SouthWest (i-1, j-1)
                if i-1 == 0 or j-1 == 0:
                    psi_SW = 0
                else:
                    psi_SW = psi[(j-1)*n + i-1]

            
                self.mv[k] = 28*psi_C - 8*(psi_N + psi_S + psi_E + psi_W) + psi_NE + psi_SE + psi_NW + psi_SW

        return self.mv

# -----------------------------------------------------------------------------
# rhs = c0*A + c1*(B - C)

# A = u_S - u_N + v_E - v_W
# B = v_C * (u_E + u_W + u_N + u_S)
# C = u_C * (v_E + v_W + v_N + v_S)

def update_rhs(tri, u, v): #
    
    n = tri.Nx 
    nc = tri.apex # i-index of triangle apex
    m = tri.Ny
    
    rhs = np.zeros(n*m)

    slope = tri.slope
    U = tri.U
    
    c0 = 3 * tri.dx
    c1 = 0.5 * tri.dx**2 * tri.Re
    
    for k in range(n*m):
        i = k % n
        j = k // n

        # boundary & exterior zones
        if i == 0 or j == 0 or i == n-1 or j == m-1:
            rhs[k] = 0
               
        elif (i < nc and j < (m-1) - slope*i) or (i > nc and j < slope * (i - nc)):
            rhs[k] = 0
            
        # interior
        else:

            # (u,v) at 9 point stencil
            #k = j*n + i
            u_C = u[k]
            v_C = v[k]
            
            # North (i, j+1)
            if j+1 == m-1: 
                u_N = U
                v_N = 0
            else:
                k_N = (j+1)*n + i
                u_N = u[k_N]
                v_N = v[k_N]
                
            # East (i+1, j)
            if i+1 == n-1: 
                u_E = 0
                v_E = 0
            else:
                k_E = j*n + i + 1
                u_E = u[k_E]
                v_E = v[k_E]
            
            # South (i, j-1)
            if j-1 == 0: 
                u_S = 0
                v_S = 0
            else:
                k_S = (j-1)*n + i
                u_S = u[k_S]
                v_S = v[k_S]

                
            # West (i-1, j)  
            if i-1 == 0: 
                u_W = 0
                v_W = 0
            else:
                k_W = j*n + i - 1
                u_W = u[k_W]
                v_W = v[k_W]
             
            A = u_S - u_N + v_E - v_W
            B = v_C * (u_E + u_W + u_N + u_S)
            C = u_C * (v_E + v_W + v_N + v_S)

            rhs[k] = c0 * A + c1 * (B - C)
            
    return rhs