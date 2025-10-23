# -*- coding: utf-8 -*-
"""
Created on Mon May 19 17:32:18 2025

@author: sarah
"""
import numpy as np

def make_rhs(height, BC): 
    n = height.N_regions-1
    rhs = np.zeros(2*n + 1)
    
    rhs[0] = -BC.p0/height.widths[0]
    rhs[n] = BC.pN/height.widths[-1]

    c = 6*BC.U #*height.visc*
    
    for k in range(n):

        rhs[n+1 + k] = (height.h_steps[k+1] - height.h_steps[k]) * c
 
    return rhs


def schurLU_solve(height, BC):
    n = height.N_regions-1
    rhs = make_rhs(height, BC)
    
    center_diag, off_diag = get_schurCompDiags(height)
    Cs, Cs_diagProd = get_Cs(n, center_diag, off_diag)
    Ds = get_Ds(n, center_diag, off_diag, Cs)

    p_peaks = np.zeros(n)
    for i in range(n):
        p_peak_ij = 0
        for j in range(n):
            # s_ij = S_ij(n, Cs_diagProd, Ds, i, j)
            if i == j:
                s_ij = Ds[i]
            elif i < j:
                s_ij = (-1)**(i+j) * Ds[i] * Cs_diagProd[i, j-1]
            else:
                s_ij = (-1)**(i+j) * Ds[j] * Cs_diagProd[j, i-1]
            
            if j == 0:
                p_peak_ij += s_ij * (rhs[n+1+j] - BC.p0 * height.h_steps[j]**3   /height.widths[j])
            elif j == n-1:
                p_peak_ij += s_ij * (rhs[n+1+j] - BC.pN * height.h_steps[j+1]**3 /height.widths[j+1] )
            else: 
                p_peak_ij += s_ij * rhs[n+1+j]
                
        p_peaks[i] = p_peak_ij
        

    p_slopes = np.zeros(n+1)
    p_slopes[0] = (p_peaks[0] - BC.p0)/height.widths[0]
    for i in range(1, n):
        p_slopes[i] = (p_peaks[i] - p_peaks[i-1])/height.widths[i]
    p_slopes[n] = (BC.pN - p_peaks[n-1])/height.widths[n]
    
    ps = make_ps(height, BC, p_slopes, p_peaks)
    return ps

def get_schurCompDiags(height):
    n = height.N_regions-1
    center_diag = np.zeros(n)
    off_diag = np.zeros(n-1)
    
    for k in range (n):
        center_diag[k] = -(height.h_steps[k]**3/height.widths[k] + height.h_steps[k+1]**3/height.widths[k+1])
        
        if k < n-1:
            off_diag[k] = height.h_steps[k+1]**3 /height.widths[k+1]
        
    return center_diag, off_diag


# recursive sequence {Si} 
def get_Cs(n, center_diag, off_diag):
    Cs = np.zeros(n-1)
    
    Cs[n-2] = off_diag[n-2] / center_diag[n-1]
    
    for k in reversed(range(n-2)):
        Cs[k] = off_diag[k] / (center_diag[k+1] - Cs[k+1]*off_diag[k+1])
    
    Cs_diagProd = np.diag(Cs)
    for k in range(n-2):
        for j in range(k+1, n-1):
            Cs_diagProd[k,j] = Cs_diagProd[k, j-1] * Cs_diagProd[j, j]
    
    return Cs, Cs_diagProd

def get_Ds(n, center_diag, off_diag, Cs):
    Ds = np.zeros(n)
    
    Ds[0] = 1/(center_diag[0] - off_diag[0]*Cs[0])
    for i in range(1, n-1):
        Ds[i] = (1 + off_diag[i-1]*Ds[i-1]*Cs[i-1]) / (center_diag[i] - off_diag[i]*Cs[i])
    Ds[n-1] = (1 + off_diag[n-2]*Ds[n-2]*Cs[n-2])/center_diag[n-1]

    return Ds

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
        
        if i >= height.i_peaks[k+1] and k < height.N_regions-1:
            k+= 1
            x_k = height.x_peaks[k]
            p_k = extrema[k-1]
            slope_k = slopes[k]
        ps[i] = slope_k*(x-x_k) + p_k

    return ps

