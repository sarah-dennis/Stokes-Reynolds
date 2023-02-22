#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 21 12:06:27 2023

@author: sarahdennis
"""
import numpy as np

n = 10
dx = 1/10
h_steps = np.random.rand(n+1)


def makeM(n, dx, h_steps):
    
    p0 = 1
    pN = 1
    
    eta, U = 1, 2

    #---------------
    M = np.zeros((2*n + 1, 2*n + 1))
    
    #B:= top left corner of M, 1/dx diagonal
    B = np.zeros((n+1, n))
    B_diag_neg = [-1/dx]*n
    B_diag_pos = [1/dx]*n
    B[0:n , 0:n] += np.diagflat(B_diag_neg)
    B[1:n+1, 0:n] += np.diagflat(B_diag_pos)  
    
    #C:= bottom right corner of M, hj-hi diagonal
    C = np.zeros((n, n+1))
    C_diag_neg = [-h**3 for h in h_steps[0:n]]
    C_diag_pos = [h**3 for h in h_steps[1:n+1]]
    C[0:n, 0:n] += np.diagflat(C_diag_neg)
    C[0:n, 1:n+1] += np.diagflat(C_diag_pos)

    #---------------
    M[0:n+1, 0:n+1] = np.identity(n+1)
    M[0:n+1, n+1:2*n+1] = B
    M[n+1:2*n+1, 0:n+1] = C
    #print("M: \n", M)
    
    #---------------

    rhs = np.zeros(2*n + 1)
    
    rhs[0] = -p0/dx
    rhs[n] = pN/dx
    
    for k in range(n):
        rhs[n+1 + k] = (h_steps[k+1] - h_steps[k]) * 6*eta*U
    
    #print("rhs: ", rhs)
    
    #---------------
    sol = np.linalg.solve(M, rhs)
    return M, rhs, sol

M, rhs, sol = makeM(n, dx, h_steps)

p_slopes = sol[0:n+1]
p_extrema =  sol[n+1:2*n+1]

