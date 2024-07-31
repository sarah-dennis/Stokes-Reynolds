#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 29 13:14:58 2024

@author: sarahdennis
"""
import numpy as np
        
from scipy.sparse.linalg import LinearOperator 
from scipy.sparse.linalg import gmres


from solvers import P_Solver


class Solver(P_Solver):   
    def __init__(self, height, p0, pN):
        p_str = "Reynolds Piecewise Analytic"
        self.rhs = make_rhs(height, p0, pN)
        self.linOp = stLinOp(height)
        super().__init__(height, p0, pN, self.solve, p_str)
    
    def solve(self, tol=1e-12):
        cs, exit_code = gmres(self.linOp, self.rhs, tol=tol)
        ps = make_ps(self.height, cs)
        return ps
    
    
def make_rhs(height, p0, pN):
    N = height.N_regions
    hs = height.h_peaks
    slopes = height.slopes
    widths = height.widths
    
    rhs = np.zeros(N+1)
    c = -6*height.visc * height.U
    
    if slopes[0] != 0:
        rhs[0] = c * (1/hs[0]) * (1/slopes[0]) - p0
    else: 
        rhs[0] = c * (-1/hs[0]**2) * widths[0] - p0
    
    for i in range(1, N):
        
        if slopes[i] != 0 and slopes[i-1] != 0:
            rhs[i] = c * (1/hs[i]) * (1/slopes[i] - 1/slopes[i-1])
            
        elif slopes[i] == 0 and slopes[i-1] != 0:
            rhs[i] = c * (-1/hs[i]**2 * widths[i] - 1/hs[i] * 1/slopes[i-1])
            
        elif slopes[i] != 0 and slopes[i-1] == 0:
            rhs[i] = c * (1/hs[i] * 1/slopes[i] -1/hs[i]**2 * widths[i-1])
            
    if slopes[N-1] != 0:   
        rhs[N] = c * (1/hs[N]) * (1/slopes[N-1]) - pN
    else:
        rhs[N] = c * (-1/hs[N]**2) * widths[N-1] - pN
        
    return rhs
        
def make_ps(height, cs):
    N = height.N_regions
    slopes = height.slopes
    widths = height.widths
    hs = height.hs
    x_peaks = height.x_peaks
    xs = height.xs
    
    ps = np.zeros(height.Nx)
    k = 0
    
    cq = -1/2 * cs[-1]
    cu = 6 * height.visc * height.U
    
    for i in range(height.Nx):
        
        if xs[i] > x_peaks[k+1] and k+1 < N:
            k += 1
        
        if slopes[k] != 0:
            dxdh = 1/slopes[k]
            hr = 1/hs[i]
            ps[i] = cq * dxdh * hr**2 - cu * dxdh * hr + cs[k]
        else:
            dx = widths[k]
            hr = 1/hs[i]
            ps[i] = cq * dx * hr**3 + cu * hr**2 * dx + cs[k]
    return ps     

              
class stLinOp(LinearOperator):
    def __init__(self, height):

        self.N = height.N_regions
        self.h_peaks = height.h_peaks
        self.widths = height.widths
        self.slopes = height.slopes
        
        self.shape = (self.N+1,self.N+1)
        self.dtype = np.dtype('f8')
        self.mv = np.zeros(self.N+1)
         
    def _matvec(self, v):
        hs = self.h_peaks
        slopes = self.slopes
        widths = self.widths
        N = self.N
        
        mv = np.zeros(N+1)
        cq = v[N]
        
        if slopes[0] != 0:
            mv[0] = -v[0] + cq * (1/(2 * hs[0]**2)) * 1/slopes[0]
        else:
            mv[0] = -v[0] + cq * 1/hs[0]**3 * widths[0]
 
        for i in range(1, N):
            if slopes[i] != 0 and slopes[i-1] != 0:
                mv[i]= v[i-1] - v[i] + cq * (1/(2 * hs[i]**2)) * (1/slopes[i] - 1/slopes[i-1])
                
            elif slopes[i] == 0 and slopes[i-1] != 0:
                mv[i] = v[i-1]-v[i] + cq * (1/hs[i]**3 * widths[i] - 1/(2 * hs[i]**2) * 1/slopes[i-1])
                
            elif slopes[i] != 0 and slopes[i-1] == 0:
                mv[i] = v[i-1]-v[i] + cq * (1/(2 * hs[i]**2) * 1/slopes[i] - 1/hs[i]**3 * widths[i-1])
    
        if slopes[N-1] != 0:
            mv[N] = -v[N-1] + cq * (1/(2 * hs[N]**2)) * 1/slopes[N-1]
        else:
            mv[N] = -v[N-1] + cq * 1/hs[N]**3 * widths[N-1]

        return mv
        
        
    
        
        
        
        
        
        
        
        
        
