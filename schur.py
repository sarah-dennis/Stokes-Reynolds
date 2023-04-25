# -*- coding: utf-8 -*-
"""
Created on Sat Apr  1 12:09:43 2023

@author: sarah
"""
import numpy as np

#------------------------------------------------------------------------------
# def schurComp(M, n, m):
#     A = M[0:n, 0:n]
#     A_inv = np.linalg.inv(A)
    
#     B1 = M[0:n, n:n+m]

#     B2 = M[n:n+m, 0:n]
    
#     C = M[n:n+m, n:n+m]

#     return -(C + np.matmul(np.matmul(B2, A_inv), B1))


#------------------------------------------------------------------------------
#build M_inv from M using schur comp
# M = [[A, B1], [B2, C]]
# M is size (n+m, n+m) 

# def build_schurComp_Minv(M, n, m):
    
#     A = M[0:n, 0:n]
#     A_inv = np.linalg.inv(A)
    
#     B1 = M[0:n, n:n+m]

#     B2 = M[n:n+m, 0:n]
    
#     C = M[n:n+m, n:n+m]

#     B2_Ainv = np.matmul(B2, A_inv)
    
#     K = -(C + np.matmul(B2_Ainv, B1))
        
#     S = np.linalg.inv(K)
    
#     M_inv = np.zeros((n+m, n+m))

#     M_inv[0:n, 0:n] = A_inv + np.matmul(A_inv, np.matmul(B1, np.matmul(S, B2_Ainv)))
#     M_inv[0:n, n:n+m] =  -np.matmul(A_inv, np.matmul(B1, S))
#     M_inv[n:n+m, 0:n] = -np.matmul(S, B2_Ainv)
#     M_inv[n:n+m, n:n+m] = S
    
#     return M_inv



def makeSymTriDiag():
    n = 3
    off_diag = [1,2]
    center_diag = [4,5,6]
    return n, off_diag, center_diag


# make schur complement K (symmetric tri-diagonal)
def make_schurCompDiags(height):
    n = height.n_steps
    hs = height.h_steps

    center_diag = np.zeros(n)
    off_diag = np.zeros(n-1)
    
    for i in range (n):
        center_diag[i] = hs[i]**3 + hs[i+1]**3
        if i < n-1:
            off_diag[i] = -hs[i+1]**3
       
    return n, off_diag, center_diag
 
#build inverse of Schur Complement using recursion
def make_SymTriInv(n, off_diag, center_diag):
   
    thetas = next_theta(np.zeros(n), 0, n, center_diag, off_diag)

    phis = next_phi(np.zeros(n), n-1, n, center_diag, off_diag)


    S = np.zeros((n,n))
    for i in range(0, n):
        for j in range(0, i+1): #j <= i
            sij = S_ij(i, j, n, thetas, phis, center_diag, off_diag)
            
            S[i,j] = sij
            
            if i != j:
                S[j, i] = sij
                
    return S


def next_theta(thetas, k, n, center_diag, off_diag):
    if k == 0:
        thetas[k] = center_diag[k]
        return next_theta(thetas, k+1, n, center_diag, off_diag)
    elif k == 1:
        thetas[k] = center_diag[k]*center_diag[k-1] - off_diag[k-1]*off_diag[k-1]
        return next_theta(thetas, k+1, n, center_diag, off_diag)
    elif k < n: 
        thetas[k] = center_diag[k]*thetas[k-1] - off_diag[k-1]*off_diag[k-1]*thetas[k-2]
        return next_theta(thetas, k+1, n, center_diag, off_diag)
    else: 
        return thetas
    
def next_phi(phis, k, n, center_diag, off_diag):
    if k == n-1:
        phis[k] = center_diag[k]
        return next_phi(phis, k-1, n, center_diag, off_diag)
    elif k == n-2:
        phis[k] = center_diag[k]* center_diag[k+1] - off_diag[k]*off_diag[k]
        return next_phi(phis, k-1, n, center_diag, off_diag)
    elif k > -1:
        phis[k] = center_diag[k]*phis[k+1] - off_diag[k]*off_diag[k]*phis[k+2]
        return next_phi(phis, k-1, n, center_diag, off_diag)
    else:
        return phis
    
    
def S_ij(i, j, n, thetas, phis, off_diag, center_diag):

    if i==j:
        if j == 0:
            return phis[i+1] / thetas[n-1]
        elif i == n-1:
            return thetas[j-1] / thetas[n-1]
        else:
            return thetas[j-1]*phis[i+1] / thetas[n-1]
    elif j < i:
        if i == n-1:
            return (-1)**(i+j) * np.prod(off_diag[j:i+1]) * thetas[j-1] / thetas[n-1]
        else:
            return (-1)**(i+j) * np.prod(off_diag[j:i+1]) * thetas[j-1]*phis[i+1] / thetas[n-1]
        