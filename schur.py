# -*- coding: utf-8 -*-
"""
Created on Wed Jun 28 10:38:12 2023

@author: sarah
"""

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

#------------------------------------------------------------------------------

#U[i,j] = prod(xs) from i to j

def triDiagProd(xs):
    U = np.diag(xs)
    n = len(xs)
    for i in range(n-1):
        for j in range(i+1, n):
            U[i,j] = U[i, j-1] * U[j, j]
    return U

# ---- RATIO METHOD -----------------------------------------------------------


# Given schur complement K (symmetric tri-diagonal)... 

# Ex. N=3
#  M = [  A1  -B1   0  ]
#      [ -B1   A2  -B2 ]
#      [  0   -B3   A3 ]

# Compute ratio elements {C_i}
# C(N-1) = B(N-1) / A(N)
# C(i) = B(i) / (A(i+1) - S(i+1) B(i+1))

def get_Cs(n, A, B):
    C = np.zeros(n-1)
    C[n-2] = B[n-2] / A[n-1]
    
    for i in reversed(range(n-2)):
        C[i] = B[i] / (A[i+1] - C[i+1]*B[i+1])
    
    return C

# Compute diagonal elements {D_i} = K^-1[i,i]

def get_Ds(n, A, B, C):
    D = np.zeros(n)
    
    D[0] = 1/(A[0] - B[0]*C[0])
    
    for i in range(1, n-1):
        D[i] = (1 + B[i-1]*D[i-1]*C[i-1]) / (A[i] - B[i]*C[i])
        
    D[n-1] = (1 + B[n-2]*D[n-2]*C[n-2])/A[n-1]

    return D
    
# Get elementwise K^-1 = S_ij
def S_ij(n, C_prod, D, i, j):
    if i == j:
        return D[i]
    # non-diagonals
    elif i < j:
        return (-1)**(i+j) * D[i] * C_prod[i, j-1]
    else:
        # j, i = i, j
        return (-1)**(i+j) * D[j] * C_prod[j, i-1]
    
def S_ij_flops(n, C, D, i, j):
    if i == j:
        return D[i]
    # non-diagonals
    elif i < j:
        return (-1)**(i+j) * D[i] * np.prod(C[i:j])
    else:
        # j, i = i, j
        return (-1)**(i+j) * D[j] * np.prod(C[j:i])
    
# Make symmetric K^-1 matrix
def get_S(n, C_prod, D):
    S = np.zeros((n,n))
    for i in range(0, n):
        for j in range(0, i+1): 
            s_ij = S_ij(n, C_prod, D, i, j)
            if i == j:
                S[i,j] = s_ij
            else:
                S[i,j] = s_ij
                S[j,i] = s_ij
            
    return S
    
    
    
    
    
    