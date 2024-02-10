# -*- coding: utf-8 -*-
"""
Created on Thu Jun  1 14:11:49 2023

@author: sarah
"""
import numpy as np
from scipy.sparse.linalg import LinearOperator

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
    
# (For testing only: Make S = K^-1 matrix)
def make_S(n, C_prod, D):
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
    
    
#------------------------------------------------------------------------------
# Reynolds equation matrix for Square Wave

# M @ lhs = rhs

# |I B1||m| = |f| <- flux continuity
# |B2 0||p|   |g| <- pressure continuity

def make_RHS(domain, height, p0, pN): 
    n = height.n_steps
    rhs = np.zeros(2*n + 1)
    
    rhs[0] = -p0/height.step_width
    rhs[n] = pN/height.step_width
    
    hs = height.h_steps
    c = 6*domain.eta*domain.U
    
    for k in range(n):
        rhs[n+1 + k] = (hs[k+1] - hs[k]) * c
    
    return rhs

# (P_extrema, P_slopes) -> [p(x)] over domain Nx
def make_ps(domain, height, p0, pN, slopes, extrema):
    ps = np.zeros(domain.Nx)
    x0 = domain.x0
    L = height.step_width

    k = 0
    x_k = domain.x0
    p_k = p0
    slope_k = slopes[k]

    for i in range(domain.Nx):
        x = domain.xs[i]
        
        #if x is in a new step
        if x > x0 + (k+1)*L:
            k += 1
            x_k = x0 + k*L
            p_k = extrema[k-1]
            slope_k = slopes[k]

        ps[i] = slope_k*(x-x_k) + p_k

    return ps

#------------------------------------------------------------------------------
# Schur Complement for square wave
#------------------------------------------------------------------------------
def make_schurCompDiags(height):
    n = height.n_steps
    hs = height.h_steps
    L = height.step_width
    
    center_diag = np.zeros(n)
    off_diag = np.zeros(n-1)
    
    for i in range (n):
        center_diag[i] = hs[i]**3 + hs[i+1]**3
        if i < n-1:
            off_diag[i] = -hs[i+1]**3
    return (-1/L) * center_diag, (-1/L) * off_diag

#------------------------------------------------------------------------------
# Helpers for pressures.squarewave_schurInvSolve()
#------------------------------------------------------------------------------

# (matrix builder only used for testing and finding cond num)
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

# Solve M_inv @ rhs = lhs
# {C_prod, D} <- Ratio Algorithm
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
            rightBlock_dotProd += S_ij(n, C_prod, D,  i - (n+1), j) * rhs[n+1 + j]
        
    return leftBlock_o + leftBlock_n + rightBlock_dotProd


#M_inv top left    
def Id_B1_schurCompInv_B2_ij (height, C_prod, D, i, j):
    L = height.step_width
    n = height.n_steps
    hj = height.h_steps[j]
    
    if i == 0 and j == 0:
        return 1 + (1/L) * hj**3 * S_ij(n, C_prod, D, i, j)
    
    elif i == 0 and j < n :
        return (-1/L) * hj**3 * (S_ij(n, C_prod, D,  i, j-1) - S_ij(n, C_prod, D, i,j))
    
    elif i == 0 and j == n: 
        return (-1/L) * hj**3 * (S_ij(n, C_prod, D,  i, j-1))
    
    elif i == n and j == 0:
        return (-1/L) * hj**3 * (S_ij(n, C_prod, D,  i-1,j))
    
    elif i == n and j < n:
         return (1/L) * hj**3 * (S_ij(n, C_prod, D, i-1, j-1) - S_ij(n, C_prod, D,  i-1, j))
    
    elif i == n and j == n:
        return 1 + (1/L) * hj**3 * S_ij(n, C_prod, D,  i-1, j-1)
    
    elif i < n and j == 0:
        return (-1/L) * hj**3 * (S_ij(n, C_prod, D,  i-1, j) - S_ij(n, C_prod, D, i, j))

    elif i < n and j == n:
        return (1/L) * hj**3 * (S_ij(n, C_prod, D,  i-1, j-1) - S_ij(n, C_prod, D,  i, j-1))
    
    else:
        return (i==j) + (1/L) * hj**3 * (S_ij(n, C_prod, D, i-1,j-1) - S_ij(n, C_prod, D,  i-1, j) - S_ij(n, C_prod, D, i, j-1) + S_ij(n, C_prod, D, i, j))    

#M_inv bottom left
def neg_schurCompInv_B2_ij(height, C_prod, D,  i, j):
    n = height.n_steps
    hj = height.h_steps[j]
    if j == 0:
        return hj**3 * S_ij(n, C_prod, D,  i,j)
    elif j < n:
        return -hj**3 * (S_ij(n, C_prod, D, i, j-1) - S_ij(n, C_prod, D,  i,j))
    else:
        return -hj**3 * S_ij(n, C_prod, D,  i, j-1)
    
#M_inv top right
def neg_B1_schurCompInv_ij(height, C_prod, D,  i, j):
    L = height.step_width
    n = height.n_steps
    if i == 0:
        return (-1/L) * (-S_ij(n,C_prod, D,  i, j))
    elif i < n:
        return (-1/L) * (S_ij(n, C_prod, D,  i-1, j) - S_ij(n, C_prod, D,  i, j))
    else:
        return (-1/L) * S_ij(n, C_prod, D, i-1, j)

#-------------------------------------------------------------
# Helpers for pressure.squarewave_pySolve()
#------------------------------------------------------------------------------
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

#------------------------------------------------------------------------------
# Helpers for pressures.squarewave_gmResSolve()
#------------------------------------------------------------------------------

class swLinOp(LinearOperator):
    #n:= number of steps 
    #L:= length of each step
    #hs:= height of each step (length n+1)
    def __init__(self, n, L, hs):
        
        self.shape = (2*n + 1, 2*n+1)
        self.L = L
        self.n = n
        self.hs = hs
        self.dtype = np.dtype('f8')
        self.mv = np.zeros(2*n+1)
        self.hscube = hs**3

        
    #M = [[I, B],[C, 0]]
    
    #[I, B][v1]
    #[C, 0][v2]
    
    #Mv = rhs

    def _matvec(self, v):
        n = self.n
        Linv = 1/self.L
        
        #----------------------
        # Upper blocks: I v1 + B v2

        #i = 0
        self.mv[0] = v[0] - v[n+1]*Linv;
        
        # 0 < i < n
        for i in range(1, n):
            self.mv[i] =  v[i] + Linv*v[i+n] + -Linv*v[i+n+1]
        # self.mv[1:n-1]=v[1:n-1] + (v[n+1:2*n-1] - v[n+2:2*n])*Linv
        
        #i = n
        self.mv[n] = v[n] + Linv*v[2*n]
        
        #---------------------- 
        #Lower blocks: C v1 + 0 v2
        # n < i < 2*n+1 
        for i in range(n+1, 2*n+1):
        
          self.mv[i] = -(self.hscube[i - n-1])*v[i-n-1] + (self.hscube[i-n])*v[i-n]
        # self.mv[n+1:2*n] = -self.hscube[0:n-1]*v[0:n-1] + self.hscube[1:n]*v[1:n]

        return self.mv