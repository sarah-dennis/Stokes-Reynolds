# -*- coding: utf-8 -*-
"""
Created on Sat Apr  1 12:09:43 2023

@author: sarah
"""
import numpy as np
# index in the notes is array indexing
# assumes as = cs

h1 = 1
h2 = 0.5
h3 = 1
h4 = 0.5
hs = [0, h1, h2, h3, h4]
n = 3

def make_S_inv(n, hs):

    #all should have length n+2 
    As = [0, -h2**3,       -h3**3,         1,             1]
    Bs = [0, h1**3 + h2**3, h2**3 + h3**3, h3**3 + h4**3, 1]
    Cs = [0, -h2**3,       -h3**3,         1,             1]
    
    thetas = next_theta(0, n, [], As, Bs, Cs)
    phis = next_phi(n+1, n, [], As, Bs, Cs)
    
    S_inv = np.zeros((n,n))
    for i in range(1, n+1):
        for j in range(i, n+1):
            
            sij = S_inv_ij(As, Bs, Cs, thetas, phis, i, j, n)
            
            S_inv[i-1,j-1] = sij
            if i != j:
                S_inv[j-1, i-1] = sij
    return S_inv
 
def next_theta(k, n, thetas, As, Bs, Cs):
    if k < 1: #k=0, k=-1
        thetas = np.zeros(n+2)
        thetas[0] = 1
        return next_theta(1, n, thetas, As, Bs, Cs)
    elif k < n+1: 
        thetas[k] = Bs[k]*thetas[k-1] - As[k]*Cs[k-1]*thetas[k-2]
        return next_theta(k+1, n, thetas, As, Bs, Cs)
    else: 
        return thetas
    
def next_phi(k, n, phis, As, Bs, Cs):
    if k > n:
        phis = np.zeros(n+3)
        phis[n+1] = 1

        return next_phi(n, n, phis, As, Bs, Cs)
    elif k > 0:
        phis[k] = Bs[k]*phis[k+1] - Cs[k]*As[k+1]*phis[k+2]
        return next_phi(k-1, n, phis, As, Bs, Cs)
    else:
        return phis
    
    
def S_inv_ij(As, Bs, Cs, thetas, phis, i, j, n):
    if i==j:
        return thetas[i-1]*phis[j+1]/thetas[n]
    elif i < j:
        return (-1)**(i+j) * np.prod(Cs[i:j]) * thetas[i-1]*phis[j+1]/thetas[n]
    else:
        return (-1)**(i+j) * np.prod(As[j:i]) * thetas[j-1]*phis[i+1]/thetas[n]
            