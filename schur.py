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
def make_S(n, off_diag, center_diag):

    thetas = get_thetas(n, off_diag, center_diag)

    phis = get_phis(n, off_diag, center_diag)
    
    off_diag_prod = triDiagProd(off_diag)
    
    S = np.zeros((n,n))
    for i in range(0, n):
        for j in range(0, i+1): #j <= i
            sij = S_ij(n, thetas, phis, off_diag_prod, i, j)
            S[i,j] = sij
            
            if i != j:
                S[j, i] = sij
    return S

    
# get (ith, jth) element of the inverse of (n, n) symmetric tri-diagonal    
#    -> excpect upper triangular i > j
def S_ij(n, thetas, phis, off_diag_prod, i, j):
    if j > i:
        j, i = i, j
        
    if i==j:
        if j == 0:
            return phis[i+1] / thetas[n-1]
        elif i == n-1:
            return thetas[j-1] / thetas[n-1]
        else:
            return thetas[j-1]*phis[i+1] / thetas[n-1]
    
    else: #(j < i)

        if j == 0 and i == n-1:
            return (-1)**(i+j) * off_diag_prod[j,i-1] / thetas[n-1]
        elif j == 0:
            return (-1)**(i+j) * off_diag_prod[j,i-1] * phis[i+1] / thetas[n-1]
        elif i == n-1:
            return (-1)**(i+j) * off_diag_prod[j,i-1] * thetas[j-1] / thetas[n-1]
        else:
            return (-1)**(i+j) * off_diag_prod[j,i-1] * thetas[j-1]*phis[i+1] / thetas[n-1]

def get_thetas(n, off_diag, center_diag):
    thetas = np.zeros(n)
    #k=0
    thetas[0] = center_diag[0]
    #k=1
    thetas[1] = center_diag[1]*center_diag[0] - off_diag[0]*off_diag[0]
    for k in range(2, n):
        
        alph = center_diag[k]*thetas[k-1]
        beta = off_diag[k-1]*off_diag[k-1]*thetas[k-2]
        thetas[k] = alph - beta
        # thetas[k] = center_diag[k]*thetas[k-1] - off_diag[k-1]*off_diag[k-1]*thetas[k-2]
    return thetas

#
def get_phis(n, off_diag, center_diag):
    phis = np.zeros(n)
    # k = n-1
    phis[n-1] = center_diag[n-1]
    # k = n-2
    phis[n-2] = center_diag[n-2]* center_diag[n-1] - off_diag[n-2]*off_diag[n-2]
    for k in range(n-3, -1, -1):
        phis[k] = center_diag[k]*phis[k+1] - off_diag[k]*off_diag[k]*phis[k+2]
    return phis


def triDiagProd(xs):
    U = np.diag(xs)
    n = len(xs)
    for i in range(n-1):
        for j in range(i+1, n):
            U[i,j] = U[i, j-1] * U[j, j]
    return U
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
    