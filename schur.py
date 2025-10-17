#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 16 12:23:38 2025

@author: sarahdennis
"""

import numpy as np
import reyn_boundary as bc

# N = height.N_regions
# height.xs = [x{0}, x{1}, ..., x{N}]
# height.widths = [x{1}-x{0}, x{2}-x{1}, ..., x{N}-x{N-1}]
# height.h_steps = [h{0}, h{1}, ..., h{N-1}] where h{k} = h(x) for x{k+1} < x < x{k}

# M : 2N-1 x 2N-1 block matix

#    M    x   =   b
# | I B ||x1| = |b1| 
# | C 0 ||x2| = |b2|

# I : N rows x N cols      identity

# B : N rows x N-1 cols    | -1/dx{0}, 0, ..., 0             |
#                          | 1/dx{1}, -1/dx{1}, 0, ..., 0    | 
#                          | 0, 1/dx{2}, -1/dx{2}, 0, ..., 0 |
#                          | ...                             |
#                          | 0, ..., 0, 1/dx{N-2}, -1/dx{N-2}|
#                          | 0, ..., 0, 1/dx{N-1}            |

# C : N-1 rows x N cols    | -h{0}^3, h{1}^3, 0, ..., 0      |
#                          | 0, -h{1}^3, h{2}^3, 0, ..., 0   |
#                          | ...                             |
#                          | 0, ..., 0, -h{N-2}^3, h{N-1}^3  |

# D : N-1 rows x N-1 cols  zeros 


#           M^-1              b   = x
# | I + B K^-1 C   -B K^-1 ||b1| = |x1|
# | -K^-1 C         K^-1   ||b2| = |x2|

# schur complement 
# K = - C B : N-1 rows x N-1 cols

def make_rhs(height, BC): 
    N = height.N_regions 
    rhs = np.zeros(2*N-1)
    
    if isinstance(BC, bc.Fixed):
        rhs[0] = -BC.p0/height.widths[0] # = dp{0} - p{1}/dx{0}
        
    elif isinstance(BC, bc.Mixed):
        h0 = height.h_steps[0]
        rhs[0] = -12*BC.Q*h0**-3 + 6*BC.U*h0**-2 # = dp{0}
        

    rhs[N-1] = BC.pN/height.widths[-1] # = dp{N-1} + p{N-1}/dx{N-1}
    
    c = 6*BC.U #*height.visc*
    
    for k in range (N-1):
        rhs[N + k] = (height.h_steps[k+1] - height.h_steps[k]) * c
 
    return rhs

def M_inv_ij(height, BC, phis, thetas, offdiag_prod, i, j):
    N = height.N_regions
    
    K_inv = K_inv_ij(height, BC, phis, thetas, offdiag_prod, i, j)
    
    if i < N and j < N: # I + B K^-1 C
        B = B_ij(height, BC, i, j)
        C = C_ij(height, i, j)
        I = 1 if i == j else 0
        return I + B * K_inv * C
    
    elif i > N and j < N: # -B K^-1
        B = B_ij(height, BC, i, j)
        return - B * K_inv
    
    elif i < N and j > N: # -K^-1 C
        C = C_ij(height, i, j)
        return - K_inv * C
    
    else: # K^-1
        return K_inv
    
    
def schur_inv_solve(height, BC):
    rhs = make_rhs(height, BC)
    phis = make_phis(height, BC)
    thetas = make_thetas(height, BC)
    offdiag_prod = make_schur_offdiag_prod(height, BC)
    
    N = height.N_regions
    xs = np.zeros(2*N - 1)
    
    for i in range(2*N-1):
        xi = 0
        for j in range(2*N-1):
            m_inv_ij = M_inv_ij(height, BC, phis, thetas, offdiag_prod, i, j)
            xi += m_inv_ij
        xs[i] = xi * rhs[i]
    
    slopes = xs[0:N]
    extrema = xs[N:]
    ps = make_ps(height, BC, slopes, extrema)
    return ps
#------------------------------------------------------------------------------

# schur complement K is size N-1 x N-1 symmetric tridiagonal

#     | b0 c0  0      ...      |
#     | a1 b1 c1  0     ...    |
#     | 0  a2 b2 c2 0     ...  |
#     | ...                    |
#     | 0 ... 0 aN-3 bN-3 cN-3 |
#     | 0   ...   0  aN-2 bN-2 |

def K_ij(height, BC, i, j): # schur complement K = - C B
    hs = height.h_steps
    ws = height.widths
    
    if i == j: #bi : center diag 
        if i == 0 and isinstance(BC, bc.Mixed):
            return -(hs[i+1]**3)/ws[i+1]
    
        else:
            return -(hs[i]**3)/ws[i] -(hs[i+1]**3)/ws[i+1]
    
    elif i+1 == j: #ci : upper diag
        return (hs[i+1]**3)/ws[i+1]
    
    elif i-1 == j: # ai : lower diag
        return (hs[i]**3)/ws[i]
    
    else:
        return 0

def K_inv_ij(height, BC, phis, thetas, offdiag_prod, i, j):    
    
    if i == j:
        return thetas[j-1] * phis[i+1] / thetas[-1]
    
    elif i > j:
        offdiag_prod_ij = offdiag_prod[j+1] / offdiag_prod[i-1] # = aj+1 * aj * ... * ai
        return (-1)**(i+j) * offdiag_prod_ij * thetas[j-1] * phis[i+1] / thetas[-1]
    
    elif i < j: # K is symmetric => K^-1 is symmetric
        return K_inv_ij(height, BC, phis, thetas, offdiag_prod, j, i)
     
# 3 components {thetas, phis, offdiag_prod} for schur complement inverse K^-1

def make_thetas(height, BC):
    N = height.N_regions
    
    thetas = np.zeros(N-1)
    
    b_0 = K_ij(height, BC, 0, 0)
    thetas[0] = b_0

    b_1 = K_ij(height, BC, 1, 1)
    a_1 = K_ij(height, BC, 1, 0) # a_1 = c_0
    thetas[1] = b_1 * b_0 - a_1**2
    
    for i in range (2, N-1):
        b_i = K_ij(height, BC, i, i)
        a_i = K_ij(height, BC, i, i-1)
        thetas[i] = b_i * thetas[i-1] -  thetas[i-2] * (a_i**2)
    return thetas
    
    
def make_phis(height, BC):
    N = height.N_regions #=n-1
    
    phis = np.zeros(N-1)
    
    b_N_minus_2 = K_ij(height, BC, N-2, N-2) #=bn
    phis[N-2] = b_N_minus_2
    
    b_N_minus_3 = K_ij(height, BC, N-3, N-3)
    c_N_minus_3 = K_ij(height, BC, N-3, N-2)
    phis[N-3] = b_N_minus_3 * b_N_minus_2 - c_N_minus_3**2
    
    for i in range(N-4, -1, -1):
        b_i = K_ij(height, BC, i, i)
        c_i = K_ij(height, BC, i, i+1)
        phis[i] = b_i * phis[i+1] - phis[i+2] * (c_i**2)
    
    return phis

def make_schur_offdiag_prod(height, BC):
    N = height.N_regions
        
    # offdiag_prod[i] = prod a1 * a2 * ... * ai 
    offdiag_prod = np.zeros(N-2)

    offdiag_prod[0] = 1
    
    for i in range(1, N-2):
        offdiag_prod[i] = K_ij(height, BC, i, i-1) * offdiag_prod[i-1]
    
    return offdiag_prod


#------------------------------------------------------------------------------

# B_ij, C_ij assume (i,j) start from (0,0) 
    
def B_ij(height, BC, i, j): # B: upper right block (N x N-1)
    ws = height.widths
    
    if i == j: #center diag
        if i == 0 and isinstance(BC, bc.Mixed):
            return 0
        else:
            return -1/ws[i]
    
    elif i-1 == j: #lower diag
        return 1/ws[i]
    
    else:
        return 0
    
def C_ij(height, i, j): # C: lower left block (N-1 x N)
    hs = height.h_steps

    if i == j: #center diag
        return -(hs[i]**3)
    
    elif i+1 == j: #upper diag
        return (hs[i+1]**3)
    
    else:
        return 0

#------------------------------------------------------------------------------


def make_ps(height, BC, slopes, extrema):
    ps = np.zeros(height.Nx)
    x0 = height.x0
    k = 0
    x_k = x0
    
    if isinstance(BC, bc.Fixed):
        p_k = BC.p0
    elif isinstance(BC, bc.Mixed):
        p_k = extrema[0] - slopes[0]*height.widths[0]
        
    slope_k = slopes[0]

    for i in range(height.Nx):
        x = height.xs[i]
        
        if i >= height.i_peaks[k+1] and k < height.N_regions-1:
            k+= 1
            x_k = height.x_peaks[k]
            p_k = extrema[k-1]
            slope_k = slopes[k]
        ps[i] = slope_k*(x-x_k) + p_k

    return ps