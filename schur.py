# -*- coding: utf-8 -*-
"""
Created on Sat Apr  1 12:09:43 2023

@author: sarah
"""
import numpy as np
#------------------------------------------------------------------------------
# Schur Complement K of (n x m block) M
#   
#  M = [ A   B1 ]
#      [ B2  C  ]

#  K = - (C + B2 A^-1 B1)

#  M^-1 = ...
#------------------------------------------------------------------------------

# build the inverse of (n, n) symmetric tri-diagonal 
def make_symTriInv(n, off_diag, center_diag):

    thetas = get_thetas(n, off_diag, center_diag)

    phis = get_phis(n, off_diag, center_diag)
    
    # off_diag_prod= get_triDiagProduct(off_diag)
    off_diag_prod= triDiagProd(off_diag)
    
    S = np.zeros((n,n))
    for i in range(0, n):
        for j in range(0, i+1): #j <= i
            sij = symTriInv_ij(i, j, n, thetas, phis, off_diag_prod)
            S[i,j] = sij
            
            if i != j:
                S[j, i] = sij
    return S

    
# get (ith, jth) element of the inverse of (n, n) symmetric tri-diagonal    
#    -> excpect upper triangular i > j
def symTriInv_ij(i, j, n, thetas, phis, off_diag_prod):
    if i==j:
        if j == 0:
            return phis[i+1] / thetas[n-1]
        elif i == n-1:
            return thetas[j-1] / thetas[n-1]
        else:
            return thetas[j-1]*phis[i+1] / thetas[n-1]
    elif j < i:
        if j == 0 and i == n-1:
            return (-1)**(i+j) * off_diag_prod[j,i-1] / thetas[n-1]
        elif j == 0:
            return (-1)**(i+j) * off_diag_prod[j,i-1] * phis[i+1] / thetas[n-1]
        elif i == n-1:
            return (-1)**(i+j) * off_diag_prod[j,i-1] * thetas[j-1] / thetas[n-1]
        else:
            return (-1)**(i+j) * off_diag_prod[j,i-1] * thetas[j-1]*phis[i+1] / thetas[n-1]

#
def get_thetas(n, off_diag, center_diag):
    return next_theta(np.zeros(n), 0, n, off_diag,  center_diag)

def next_theta(thetas, k, n, off_diag, center_diag):
    if k == 0:
        thetas[k] = center_diag[k]
        return next_theta(thetas, k+1, n, off_diag, center_diag)
    elif k == 1:
        thetas[k] = center_diag[k]*center_diag[k-1] - off_diag[k-1]*off_diag[k-1]
        return next_theta(thetas, k+1, n, off_diag, center_diag)
    elif k < n: 
        thetas[k] = center_diag[k]*thetas[k-1] - off_diag[k-1]*off_diag[k-1]*thetas[k-2]
        return next_theta(thetas, k+1, n, off_diag, center_diag)
    else: 
        return thetas

# 
def get_phis(n, off_diag, center_diag):
    return next_phi(np.zeros(n), n-1, n, off_diag, center_diag)
    
def next_phi(phis, k, n, off_diag, center_diag):
    if k == n-1:
        phis[k] = center_diag[k]
        return next_phi(phis, k-1, n, off_diag, center_diag)
    elif k == n-2:
        phis[k] = center_diag[k]* center_diag[k+1] - off_diag[k]*off_diag[k]
        return next_phi(phis, k-1, n, off_diag, center_diag)
    elif k > -1:
        phis[k] = center_diag[k]*phis[k+1] - off_diag[k]*off_diag[k]*phis[k+2]
        return next_phi(phis, k-1, n, off_diag, center_diag)
    else:

        return phis
    
def triDiagProd(xs):
    U = np.diag(xs)
    n = len(xs)
    for i in range(n-1):
        for j in range(i+1, n):
            U[i,j] = U[i, j-1] * U[j, j]
         
    return U
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
    