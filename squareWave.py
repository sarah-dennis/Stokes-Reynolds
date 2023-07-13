# -*- coding: utf-8 -*-
"""
Created on Thu Jun  1 14:11:49 2023

@author: sarah
"""
import numpy as np

import schur

# M @ lhs = rhs

# |I B1||m| = |f|
# |B2 0||p|   |g|

def make_RHS(domain, height, p0, pN): 
    n = height.n_steps
    rhs = np.zeros(2*n + 1)
    
    rhs[0] = -p0/height.step_width
    rhs[n] = pN/height.step_width
    
    for k in range(n):
        rhs[n+1 + k] = (height.h_steps[k+1] - height.h_steps[k]) * 6*domain.eta*domain.U
        
    return rhs

# Make diagonals of schur complement K
def make_schurCompDiags(height):
    n = height.n_steps
    hs = height.h_steps

    center_diag = np.zeros(n)
    off_diag = np.zeros(n-1)
    
    for i in range (n):
        center_diag[i] = hs[i]**3 + hs[i+1]**3
        if i < n-1:
            off_diag[i] = -hs[i+1]**3
    return (-1/height.step_width) * center_diag, (-1/height.step_width) * off_diag

#------------------------------------------------------------------------------
# Helpers for pressures.schur_InvSolve()

# Solve M_inv @ rhs = lhs
def schurInvSol_i(rhs, height, C_prod, D, i):
    n = height.n_steps
    
    if i < n+1:       # solving for M_i

        leftBlock_o = Id_B1_schurCompInv_B2_ij(height, C_prod, D,  i, 0) * rhs[0]
        leftBlock_n = Id_B1_schurCompInv_B2_ij(height, C_prod, D,  i, n) * rhs[n]
        
        rightBlock_dotProd = 0
        for j in range(n):
            rightBlock_dotProd += neg_B1_schurCompInv_ij(height, C_prod, D, i, j) * rhs[n+1+j]
        
    else: #n < i < 2*n+1  #solving for P_i  

        leftBlock_o = neg_schurCompInv_B2_ij(height, C_prod, D,  i - (n+1), 0) * rhs[0]
        leftBlock_n = neg_schurCompInv_B2_ij(height, C_prod, D,  i - (n+1), n) * rhs[n]
        
        rightBlock_dotProd = 0
        for j in range(n):
            rightBlock_dotProd += schur.S_ij(n, C_prod, D,  i - (n+1), j) * rhs[n+1 + j]
        
    return leftBlock_o + leftBlock_n + rightBlock_dotProd


#M_inv top left    
def Id_B1_schurCompInv_B2_ij (height, C_prod, D, i, j):
    L = height.step_width
    n = height.n_steps
    hj = height.h_steps[j]
    
    if i == 0 and j == 0:
        return 1 + (1/L) * hj**3 * schur.S_ij(n, C_prod, D, i, j)
    
    elif i == 0 and j < n :
        return (-1/L) * hj**3 * (schur.S_ij(n, C_prod, D,  i, j-1) - schur.S_ij(n, C_prod, D, i,j))
    
    elif i == 0 and j == n: 
        return (-1/L) * hj**3 * (schur.S_ij(n, C_prod, D,  i, j-1))
    
    elif i == n and j == 0:
        return (-1/L) * hj**3 * (schur.S_ij(n, C_prod, D,  i-1,j))
    
    elif i == n and j < n:
         return (1/L) * hj**3 * (schur.S_ij(n, C_prod, D, i-1, j-1) - schur.S_ij(n, C_prod, D,  i-1, j))
    
    elif i == n and j == n:
        return 1 + (1/L) * hj**3 * schur.S_ij(n, C_prod, D,  i-1, j-1)
    
    elif i < n and j == 0:
        return (-1/L) * hj**3 * (schur.S_ij(n, C_prod, D,  i-1, j) - schur.S_ij(n, C_prod, D, i, j))

    elif i < n and j == n:
        return (1/L) * hj**3 * (schur.S_ij(n, C_prod, D,  i-1, j-1) - schur.S_ij(n, C_prod, D,  i, j-1))
    
    else:
        return (i==j) + (1/L) * hj**3 * (schur.S_ij(n, C_prod, D, i-1,j-1) - schur.S_ij(n, C_prod, D,  i-1, j) - schur.S_ij(n, C_prod, D, i, j-1) + schur.S_ij(n, C_prod, D, i, j))    

#M_inv bottom left
def neg_schurCompInv_B2_ij(height, C_prod, D,  i, j):
    n = height.n_steps
    hj = height.h_steps[j]
    if j == 0:
        return hj**3 * schur.S_ij(n, C_prod, D,  i,j)
    elif j < n:
        return -hj**3 * (schur.S_ij(n, C_prod, D, i, j-1) - schur.S_ij(n, C_prod, D,  i,j))
    else:
        return -hj**3 * schur.S_ij(n, C_prod, D,  i, j-1)
    
#M_inv top right
def neg_B1_schurCompInv_ij(height, C_prod, D,  i, j):
    L = height.step_width
    n = height.n_steps
    if i == 0:
        return (-1/L) * (-schur.S_ij(n,C_prod, D,  i, j))
    elif i < n:
        return (-1/L) * (schur.S_ij(n, C_prod, D,  i-1, j) - schur.S_ij(n, C_prod, D,  i, j))
    else:
        return (-1/L) * schur.S_ij(n, C_prod, D, i-1, j)

#------- Matrix Builders ------------------------------------------------------

def make_M(domain, height, p0, pN):

    n = height.n_steps
    hs = height.h_steps 
    #---------------
    #M = [[I, B], [C, 0]]
    M = np.zeros((2*n + 1, 2*n + 1))


    #B:= top right corner of M, 1/L diagonal
    B = np.zeros((n+1, n))
    B_diag_neg = [-1/height.step_width]*n
    B_diag_pos = [1/height.step_width]*n
    B[0:n , 0:n] += np.diagflat(B_diag_neg)
    B[1:n+1, 0:n] += np.diagflat(B_diag_pos)  
    M[0:n+1, n+1:2*n+1] = B
    
    #C:= bottom left corner of M, hj-hi diagonal
    C = np.zeros((n, n+1))
    C_diag_neg = [-h**3 for h in hs[0:n]] 
    C_diag_pos = [h**3 for h in hs[1:n+1]]
    C[0:n, 0:n] += np.diagflat(C_diag_neg)
    C[0:n, 1:n+1] += np.diagflat(C_diag_pos)
    
    #---------------
    M[0:n+1, 0:n+1] = np.identity(n+1)
    M[0:n+1, n+1:2*n+1] = B
    M[n+1:2*n+1, 0:n+1] = C
        
    return M

def make_Minv_schurComp(height, S):
    
    n = height.n_steps
    
    #Minv = [[A, B],[C, S]]
    M_inv = np.zeros((2*n + 1, 2*n + 1))
    
    A = np.zeros((n+1, n+1))
    B = np.zeros((n+1, n))
    C = np.zeros((n, n+1))
    
    for i in range(0, n+1):
        for j in range(0, n+1):
            
            A[i,j] = Id_B1_schurCompInv_B2_ij(height, S, i, j)
            
            if i < n:
                C[i, j] = neg_schurCompInv_B2_ij(height, S, i, j)
            
            if j < n:
                B[i,j] = neg_B1_schurCompInv_ij(height, S, i, j)
                
    M_inv[0:n+1, 0:n+1] = A
    M_inv[0:n+1, n+1:2*n+1] = B
    M_inv[n+1:2*n+1, 0:n+1] = C
    M_inv[n+1:2*n+1, n+1:2*n+1] = S
            
    return M_inv

# -----------------------------------------------------------------------------
# Construct piecewise linear pressure on Nx grid from list of extrema
def make_ps(domain, height, p0, pN, slopes, extrema):
    ps = np.zeros(domain.Nx)
    L = height.step_width
    
    k = 0
    x_k = domain.x0
    p_k = p0
    slope_k = slopes[k]

    for i in range(domain.Nx-1):
        x = domain.xs[i]
        
        #if x is in a new step
        if x > domain.x0 + (k+1)*L:
            k += 1
            x_k = domain.x0 + k*L
            p_k = extrema[k-1]
            slope_k = slopes[k]

        ps[i] = slope_k*(x-x_k) + p_k
    
    ps[-1] = pN
    return ps