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

# size 2N-1 block matix

# | I B ||x1| = |b1| 
# | C 0 ||x2| = |b2|

# I : N rows x N cols     identity

# B : N rows x N-1 cols   [[-1/dx{0}, 0, ..., 0], 
#                          [1/dx{1}, -1/dx{1}, 0, ..., 0], 
#                          [0, 1/dx{2}, -1/dx{2}, 0, ..., 0]
#                       ...[0, ..., 0, 1/dx{N-2}, -1/dx{N-2}]
#                          [0, ..., 0, 1/dx{N-1}]

# C : N-1 rows x N cols   [[-h{0}^3, h{1}^3, 0, ..., 0],
#                          [0, -h{1}^3, h{2}^3, 0, ..., 0],
#                       ...[0, ..., 0, -h{N-2}^3, h{N-1}^3]

# D : N-1 rows x N-1 cols  Zeros 


# L U decomposition

# for Mx = b
# 1. Ly = b (solve for y using fwd sub)
# 2. Ux = y (solve for x using bck sub)

#    L      U    x      b
# | I 0 || I B ||x1| = |b1| 
# | C K || 0 I ||x2| = |b2|



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

def LU_solve(height, BC):
    rhs = make_rhs(height, BC)
    y = fwd_sub(height, BC, rhs)
    x = bck_sub(height, BC, y)
    
    N = height.N_regions
    p_slopes = x[0 : N]
    p_extrema = x[N : 2*N-1]
    
    ps = make_ps(height, BC, p_slopes, p_extrema)
    return ps


def fwd_sub(height, BC, rhs):
    N = height.N_regions
    
    y = np.zeros(2*N -1)
    
    for i in range (0, 2*N-1, 1):
        
        z = rhs[i]
        for j in range (i):
            z -= L_ij(height, BC, i, j) * y[j]
            
        y[i] = z / L_ij(height, BC, i, i)

    return y

def bck_sub(height, BC, y):
    N = height.N_regions
    x = np.zeros(2*N -1)
    
    for i in range(2*N-2, -1, -1):
        
        z = y[i]
        for j in range(i+1, 2*N-1):
            z -= U_ij(height, BC, i, j) * x[j]
            
        x[i] = z/U_ij(height, BC, i, i)
        
    return x

def L_ij (height, BC, i, j):
    N = height.N_regions
    
    if i < N: # upper | I  0 |
        if i == j:
            return 1
        else:
            return 0
        
    else:  # lower | C K |
        if j < N:
            return C_ij(height, i-N, j)
        else:
            return K_ij(height, BC, i-N, j)

def U_ij(height, BC, i, j):
    N = height.N_regions
    
    if i < N: # upper | I B |
        if j < N: 
            if i == j:
                return 1
            else:
                return 0
        else:
            return B_ij(height, BC, i, j-N)
        
    else: #lower | 0 I |
        if i == j:
            return 1
        else:
            return 0
        

# K_ij, B_ij, C_ij take **shifted indices** starting from 0!! 

def K_ij(height, BC, i, j): # schur complement K = - C B
    hs = height.h_steps
    ws = height.widths
    
    if i == j: #center diag
        if i == 0 and isinstance(BC, bc.Mixed):
            return -(hs[i+1]**3)/ws[i+1]
        
        else:
            return -(hs[i]**3)/ws[i] -(hs[i+1]**3)/ws[i+1]
    
    elif i+1 == j: #upper diag
        return (hs[i+1]**3)/ws[i+1]
    
    elif i-1 == j: #lower diag
        return (hs[i]**3)/ws[i]
    
    else:
        return 0
    
    
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