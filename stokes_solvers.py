# -*- coding: utf-8 -*-
"""
Created on Tue May 21 13:58:57 2024

@author: sarah
"""
import numpy as np
from scipy.sparse.linalg import LinearOperator
from scipy.sparse import csc_matrix

# -----------------------------------------------------------------------------
# rhs = c0*A + c1*(B - C)

# A = u_S - u_N + v_E - v_W
# B = v_C * (u_E + u_W + u_N + u_S)
# C = u_C * (v_E + v_W + v_N + v_S)

def update_rhs(tri, u, v): #
    
    n = tri.Nx 
    nc = n//2 # i-index of triangle apex
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

# -----------------------------------------------------------------------------

# Velocity <-> Stream update

# u[k] = c2 * (psi_N - psi_S) - c3 * (u_N + u_S)
# v[k] = -c2 * (psi_E - psi_W) - c3 * (v_E + v_W)

def uv_approx(tri, u, v, psi):
    
    n = tri.Nx
    m = tri.Ny

    nc = n//2
    slope = tri.slope
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
        elif j == 0 or i == 0 or i == n-1:
            u[k] = 0
            v[k] = 0 
        
        elif (i < nc and j < (m-1) - slope*i) or (i > nc and j < slope * (i - nc)):
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
                psi_N = psi[k_N]
                
            # East (i+1, j)
            if i+1 == n-1 : 
                v_E = 0
                psi_E = 0 
            else:
                k_E = j*n + i + 1
                v_E = v[k_E]
                psi_E = psi[k_E] 
            
            # South (i, j-1)
            if j-1 == 0 : 
                u_S = 0
                psi_S = 0 
            else:
                k_S = (j-1)*n + i
                u_S = u[k_S]
                psi_S = psi[k_S]
                
            # West (i-1, j)  
            if i-1 == 0 :
                v_W = 0
                psi_W = 0 
            else:
                k_W = j*n + i - 1
                v_W = v[k_W]
                psi_W = psi[k_W]
   
            u[k] = c2 * (psi_N - psi_S) - c3 * (u_N + u_S)
            v[k] = -c2 * (psi_E - psi_W) - c3 * (v_E + v_W)
    
    return u, v

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

def Dpsi_cscmatrixBuild(tri):
    m = tri.Ny
    n = tri.Nx

    slope = tri.slope
    nc = n//2

    data = np.zeros(m*n*9)
    row = np.zeros(m*n*9)
    col = np.zeros(m*n*9)
    
    s=0
    
    for k in range(m*n):
    
        i = k % n
        j = k // n
        
        # k = j*n + i
        
        if i == 0 or j == 0 or i == n-1 or j == m-1:
            row[s] = k
            col[s] = k
            data[s] = 1
            s+=1
        
        elif  (i < nc and j < (m-1) - slope*i) or (i > nc and j < slope * (i - nc)):
            row[s] = k
            col[s] = k
            data[s] = 1
            s+=1
        else: 
            row[s] = k
            col[s] = k
            data[s] = 28
            s+=1
            
            row[s] = k
            col[s] = j*n + i - 1
            data[s] = -8
            s+=1 
            
            row[s] = k
            col[s] = j*n + i + 1
            data[s] = -8
            s+=1 

            row[s] = k
            col[s] = (j-1)*n + i
            data[s] = -8
            s+=1
            
            row[s] = k
            col[s] = (j+1)*n + i
            data[s] = -8
            s+=1
            
            row[s] = k
            col[s] = (j-1)*n + i-1
            data[s] = 1
            s+=1
            
            row[s] = k
            col[s] = (j-1)*n + i+1
            data[s] = 1
            s+=1
            
            row[s] = k
            col[s] = (j+1)*n + i-1
            data[s] = 1
            s+=1
            
            row[s] = k
            col[s] = (j+1)*n + i+1
            data[s] = 1
            s+=1
            
    csc_mat = csc_matrix((data, (row, col)), (m*n, m*n))
    return csc_mat


# -----------------------------------------------------------------------------
# boundary mirroring

def mirror_boundary(tri, psi):
    n = tri.Nx 
    m = tri.Ny 
    slope = tri.slope
    dx = tri.dx
    
    i_mid = n//2
    
    for j in range(m-1):
        
        dj = j % slope
        di = int(j//slope)
        
        k_left = j*n + i_mid - di
        k_right = j*n + i_mid + di
        
        if dj == 0: #true boundary point
            psi[k_left] = 0 
            psi[k_right] = 0
        
        else:
    
            psi[k_left-1] = (1 - dj*dx/slope)*psi[k_left]
            psi[k_right+1] = (1 + dj*dx/slope)*psi[k_right] 
             
    return psi
        
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

