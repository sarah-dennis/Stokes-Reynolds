# -*- coding: utf-8 -*-
"""
Created on Thu Jun  1 14:11:49 2023

@author: sarah
"""
import numpy as np

from solvers import P_Solver



class Solver_numpy(P_Solver):
    def __init__(self, height, p0, pN):
        p_str = "Reynolds Piecewise Analytic (numpy.linalg)"
        self.N_steps = height.N_steps
        self.rhs = make_RHS(height, p0, pN)
        self.swMat = make_M(height, p0, pN)  
        
        super().__init__(height, p0, pN, self.numpy_solve, p_str)
    
    def numpy_solve(self):
        sol = np.linalg.solve(self.swMat, self.rhs)
        p_slopes = sol[0:self.N_steps+1]
        p_extrema = sol[self.N_steps+1:2*self.N_steps+1]
        ps = make_ps(self.height, self.p0, self.pN, p_slopes, p_extrema)
        return ps

class Solver_schurInv(P_Solver):
    def __init__(self, height, p0, pN):
        p_str = "Reynolds Piecewise Analytic (schur inv)"
        self.rhs = make_RHS(height, p0, pN)
        
        center_diag, off_diag = make_schurCompDiags(height)
        C = get_Cs(height.N_steps, center_diag, off_diag)
        self.C_prod = triDiagProd(C)
        self.D = get_Ds(height.N_steps, center_diag, off_diag, C)
        super().__init__(height, p0, pN, self.schurInv_solve, p_str)
        
    def schurInv_solve(self):
        n = self.height.N_steps

        sol = np.zeros(2*n+1)
        for i in range(2*n+1):
            sol[i] = schurInvSol_i(self.height, self.rhs, self.C_prod, self.D,  i)
        p_slopes = sol[0:n+1]
        p_extrema = sol[n+1:2*n+1]
        
        ps = make_ps(self.height, self.p0, self.pN, p_slopes, p_extrema) 
        
        return ps



class Solver_schurLU(P_Solver):
    def __init__(self, height, p0, pN):
        p_str = "Reynolds Piecewise Analytic (schur-LU)"

        self.rhs = make_RHS(height, p0, pN)
        
        center_diag, off_diag = make_schurCompDiags(height)
        C = get_Cs(height.N_steps, center_diag, off_diag)
        self.C_prod = triDiagProd(C)
        self.D = get_Ds(height.N_steps, center_diag, off_diag, C)
        
        super().__init__(height, p0, pN, self.schurLU_solve, p_str)
        
    def schurLU_solve(self):
        n = self.height.N_steps
        L = self.height.step_width
        hs = self.height.h_steps

        D = self.D
        C_prod = self.C_prod
        rhs = self.rhs
        p0 = self.p0
        pN = self.pN
        # L block -  fwd sub

        
        w = np.zeros(n)
        for i in range(n):
            w_ij = 0
            for j in range(n):
                s_ij = S_ij(n, C_prod, D, i, j)
                
                if j == 0:
                    w_ij += s_ij * (rhs[n+1+j] - (1/L) * p0 * hs[j]**3)
                elif j == n-1:
                    w_ij += s_ij * (rhs[n+1+j] -  (1/L) * pN * hs[j+1]**3)
                else: 
                    w_ij += s_ij * rhs[n+1+j]
            w[i] = w_ij
            
        # U block  - back sub

        x = np.zeros(n+1)
        
        x[0] = 1/L * (w[0] - p0)
        for i in range(1, n):
            x[i] = 1/L * (w[i] - w[i-1])
            
        x[n] = 1/L * (pN - w[n-1])
        
        ps = make_ps(self.height, p0, pN, x, w)
        return ps


#------------------------------------------------------------------------------
# Reynolds equation matrix for Square Wave

def make_RHS(height, p0, pN): 
    n = height.N_steps
    rhs = np.zeros(2*n + 1)
    
    rhs[0] = -p0/height.step_width
    rhs[n] = pN/height.step_width
    
    hs = height.h_steps
    c = 6*height.visc*height.U
    
    for k in range(n):
        rhs[n+1 + k] = (hs[k+1] - hs[k]) * c
    
    return rhs

# (P_extrema, P_slopes) -> [p(x)] over domain Nx
def make_ps(height, p0, pN, slopes, extrema):
    ps = np.zeros(height.Nx)
    x0 = height.x0
    L = height.step_width

    k = 0
    x_k = x0
    p_k = p0
    slope_k = slopes[k]

    for i in range(height.Nx):
        x = height.xs[i]
        
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

# Schur Complement K of (n x m block) M
#   
#  M = [ A   B1 ]
#      [ B2  C  ]

#  K = - (C + B2 A^-1 B1)

def make_schurCompDiags(height):
    n = height.N_steps
    hs = height.h_steps
    L = height.step_width
    
    center_diag = np.zeros(n)
    off_diag = np.zeros(n-1)
    
    for i in range (n):
        center_diag[i] = hs[i]**3 + hs[i+1]**3
        if i < n-1:
            off_diag[i] = -hs[i+1]**3
    return (-1/L) * center_diag, (-1/L) * off_diag

# -----------------------------------------------------------------------------
# Ratio Method : Tridiagonal Schur Complement -> Inverse and LU decompositions
# -----------------------------------------------------------------------------

#  K = [  A1  -B1   0  ]
#      [ -B1   A2  -B2 ]
#      [  0   -B3   A3 ]

# Compute ratio elements {C_i}
#   C(N-1) = B(N-1) / A(N)
#   C(i) = B(i) / (A(i+1) - S(i+1) B(i+1))

def get_Cs(n, A, B):
    C = np.zeros(n-1)
    C[n-2] = B[n-2] / A[n-1]
    
    for i in reversed(range(n-2)):
        C[i] = B[i] / (A[i+1] - C[i+1]*B[i+1])
    
    return C
    
def triDiagProd(C):
    diagProd = np.diag(C)
    n = len(C)
    for i in range(n-1):
        for j in range(i+1, n):
            diagProd[i,j] = diagProd[i, j-1] * diagProd[j, j]
    return diagProd


# Compute diagonal elements {D_i}
#   D(i) = K^-1[i,i]
def get_Ds(n, A, B, C):
    D = np.zeros(n)
    
    D[0] = 1/(A[0] - B[0]*C[0])
    
    for i in range(1, n-1):
        D[i] = (1 + B[i-1]*D[i-1]*C[i-1]) / (A[i] - B[i]*C[i])
        
    D[n-1] = (1 + B[n-2]*D[n-2]*C[n-2])/A[n-1]

    return D
    
# Elementwise accessing  K^-1 = S_ij
def S_ij(n, C_prod, D, i, j):
    if i == j:
        return D[i]
    # non-diagonals
    elif i < j:
        return (-1)**(i+j) * D[i] * C_prod[i, j-1]
    else:
        # j, i = i, j
        return (-1)**(i+j) * D[j] * C_prod[j, i-1]
    


#------------------------------------------------------------------------------
# Helpers for pressures.stepwave_schurInvSolve()
#------------------------------------------------------------------------------

# Solve M_inv @ rhs = lhs

def schurInvSol_i(height, rhs, C_prod, D, i):
    n = height.N_steps
    
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
    n = height.N_steps
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
    n = height.N_steps
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
    n = height.N_steps
    if i == 0:
        return (-1/L) * (-S_ij(n,C_prod, D,  i, j))
    elif i < n:
        return (-1/L) * (S_ij(n, C_prod, D,  i-1, j) - S_ij(n, C_prod, D,  i, j))
    else:
        return (-1/L) * S_ij(n, C_prod, D, i-1, j)

#------------------------------------------------------------------------------
# Helpers for pressure.squarewave_pySolve()
#------------------------------------------------------------------------------

def make_M(height, p0, pN):

    n = height.N_steps
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
