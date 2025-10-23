# -*- coding: utf-8 -*-
"""
Created on Mon May 19 17:32:18 2025

@author: sarah
"""
import numpy as np
import reyn_boundary as bc

# def make_rhs(height, BC): 
#     n = height.N_regions-1
#     rhs = np.zeros(2*n + 1)
    
#     rhs[0] = -BC.p0/height.widths[0]
#     rhs[n] = BC.pN/height.widths[-1]

#     c = 6*BC.U #*height.visc*
    
#     for k in range(n):

#         rhs[n+1 + k] = (height.h_steps[k+1] - height.h_steps[k]) * c
 
#     return rhs

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

def schurLU_solve(height, BC):
    n = height.N_regions-1
    rhs = make_rhs(height, BC)
    
    S, S_prod = get_S(height, BC)
    D = get_D(height, BC, S)
    
    p_peaks = np.zeros(n)
    for i in range(n):
        p_peak_ij = 0
        
        for j in range(n):
            if i == j:
                k_inv_ij = D[i]
                
            elif i < j:
                k_inv_ij = D[i] * (S_prod[j-1]/S_prod[i-1]) 
                
            else:
                k_inv_ij = D[j] * (S_prod[i-1]/S_prod[j-1])
            
            if j == 0:
                p_peak_ij += k_inv_ij * (rhs[n+1+j] - BC.p0 * height.h_steps[j]**3   /height.widths[j])
            elif j == n-1:
                p_peak_ij += k_inv_ij * (rhs[n+1+j] - BC.pN * height.h_steps[j+1]**3 /height.widths[j+1])
            else: 
                p_peak_ij += k_inv_ij * rhs[n+1+j]
                
        p_peaks[i] = p_peak_ij
        

    p_slopes = np.zeros(n+1)
    p_slopes[0] = (p_peaks[0] - BC.p0)/height.widths[0]
    for i in range(1, n):
        p_slopes[i] = (p_peaks[i] - p_peaks[i-1])/height.widths[i]
    p_slopes[n] = (BC.pN - p_peaks[n-1])/height.widths[n]
    
    ps = make_ps(height, BC, p_slopes, p_peaks)
    return ps



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


# recursive sequence {Si} 
def get_S(height, BC):
    N = height.N_regions
    n = N-1 # schur complement is size N-1 x N-1
    
    S = np.zeros(n-1)


    k=n-2
    off_diag = K_ij(height, BC, k, k+1)
    center_diag = K_ij(height, BC, k+1, k+1)
    # S[k] = off_diag / center_diag
    S[k] = -off_diag / center_diag
    
    for k in range(n-3, -1, -1):
        off_diag_succ = off_diag #=K_ij(height, BC, k+1, k+2)
        off_diag = K_ij(height, BC, k, k+1)
        center_diag = K_ij(height, BC, k+1, k+1)
        
        S[k] = -off_diag / (center_diag + S[k+1] * off_diag_succ)
        
        
    S_prod = np.zeros(n)
    S_prod[0] = S[0]
    for k in range(1, n-1):
        S_prod[k] = S_prod[k-1]*S[k]
    S_prod[n-1]=1

    return S, S_prod

# recursive sequence {Di} = diags of schur inverse
def get_D(height, BC, S):
    N = height.N_regions
    n = N-1 # schur complement is size N-1 x N-1
    
    D = np.zeros(n)
    
    k=0
    off_diag = K_ij(height, BC, k, k+1)
    center_diag = K_ij(height, BC, k, k)
    D[k] = 1/(center_diag + off_diag * S[k])
    
    
    for k in range(1, n-1):
        off_diag_pred = off_diag
        off_diag = K_ij(height, BC, k, k+1)
        center_diag = K_ij(height, BC, k, k)
        
        D[k] = (1 - off_diag_pred*D[k-1]*S[k-1]) / (center_diag + off_diag * S[k])
        

    off_diag_pred = off_diag
    center_diag = K_ij(height, BC, n-1, n-1)
    D[n-1] = (1 - off_diag_pred*D[n-2]*S[n-2]) / center_diag 
    
    return D

# (P_extrema, P_slopes) -> [p(x)] over domain Nx
def make_ps(height, BC, slopes, extrema):
    ps = np.zeros(height.Nx)
    x0 = height.x0
    k = 0
    x_k = x0
    p_k = BC.p0
    slope_k = slopes[0]

    for i in range(height.Nx):
        x = height.xs[i]
        ps[i] = slope_k*(x-x_k) + p_k
        
        if i >= height.i_peaks[k+1] and k < height.N_regions-1:
            k+= 1
            x_k = height.x_peaks[k]
            p_k = extrema[k-1]
            slope_k = slopes[k]
        

    return ps

