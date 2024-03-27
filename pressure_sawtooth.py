#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 29 13:14:58 2024

@author: sarahdennis
"""

        
from scipy.sparse.linalg import LinearOperator 
import numpy as np

class stLinOp(LinearOperator):
    def __init__(self, N, hs, slopes):
        self.shape = (N+1,N+1)
        self.dtype = np.dtype('f8')
        self.N = N
        self.hs = hs
        
        self.slopes = slopes

        self.mv = np.zeros(N+1)

        
    def _matvec(self, v):
        hs = self.hs
        slopes = self.slopes
        N = self.N
        
        mv = np.zeros(N+1)
        cq = v[N]
        
        mv[0] = -v[0] + cq * (1/(2 * hs[0]**2)) * 1/slopes[0]

        for i in range(1, N):
            mv[i]= v[i-1] - v[i] + cq * (1/(2 * hs[i]**2)) * (1/slopes[i] - 1/slopes[i-1])
        mv[N] = -v[N-1] + cq * (1/(2 * hs[N]**2)) * 1/slopes[N-1]

        return mv
            
            
def make_rhs(N, hs, slopes, p0, pN, eta_U):
    rhs = np.zeros(N+1)
    c = -6*eta_U
    rhs[0] = c * (1/hs[0]) * (1/slopes[0]) - p0
    
    for i in range(1, N):
        rhs[i] = c * (1/hs[i]) * (1/slopes[i] - 1/slopes[i-1])
    rhs[N] = c * (1/hs[N]) * (1/slopes[N-1]) - pN
    return rhs

def make_ps(domain, height, vs, N):
    ps = np.zeros(domain.Nx)
    slopes = height.slopes
    hs = height.hs
    x_peaks = height.x_peaks
    xs = domain.xs
    
    k = 0
    
    cq = -1/2 * vs[-1]
    cu = -6 * domain.eta * domain.U
    
    for i in range(domain.Nx):
        
        if xs[i] > x_peaks[k+1] and k+1 < N:
            k += 1
        
        dxdh = 1/slopes[k]
        hr = 1/hs[i]
        ps[i] = cq * dxdh * hr**2 + cu * dxdh * hr + vs[k]
    return ps
        
        
        
class stLinOp_two(LinearOperator):
    def __init__(self, N, hs, slopes, dxs):
        self.shape = (N+1,N+1)
        self.dtype = np.dtype('f8')
        self.N = N
        self.hs = hs
        self.dxs = dxs
         
        self.slopes = slopes
    
        self.mv = np.zeros(N+1)
         
    def _matvec(self, v):
        hs = self.hs
        slopes = self.slopes
        dxs = self.dxs
        N = self.N
        
        mv = np.zeros(N+1)
        cq = v[N]
        
        if slopes[0] != 0:
            mv[0] = -v[0] + cq * (1/(2 * hs[0]**2)) * 1/slopes[0]
        else:
            mv[0] = -v[0] + cq * 1/hs[0]**3 * dxs[0]
 
        for i in range(1, N):
            if slopes[i] != 0 and slopes[i-1] != 0:
                mv[i]= v[i-1] - v[i] + cq * (1/(2 * hs[i]**2)) * (1/slopes[i] - 1/slopes[i-1])
                
            elif slopes[i] == 0 and slopes[i-1] != 0:
                mv[i] = v[i-1]-v[i] + cq * (1/hs[i]**3 * dxs[i] - 1/(2 * hs[i]**2) * 1/slopes[i-1])
                
            elif slopes[i] != 0 and slopes[i-1] == 0:
                mv[i] = v[i-1]-v[i] + cq * (1/(2 * hs[i]**2) * 1/slopes[i] - 1/hs[i]**3 * dxs[i])
    
        if slopes[N-1] != 0:
            mv[N] = -v[N-1] + cq * (1/(2 * hs[N]**2)) * 1/slopes[N-1]
        else:
            mv[N] = -v[N-1] + cq * 1/hs[N]**3 * dxs[N-1]

        return mv
        
        
    
def make_rhs_two(N, hs, slopes, dxs, p0, pN, eta_U):
    rhs = np.zeros(N+1)
    c = -6*eta_U
    
    if slopes[0] != 0:
        rhs[0] = c * (1/hs[0]) * (1/slopes[0]) - p0
    else: 
        rhs[0] = c * (-1/hs[0]**2) * dxs[0] - p0
    
    for i in range(1, N):
        
        if slopes[i] != 0 and slopes[i-1] != 0:
            rhs[i] = c * (1/hs[i]) * (1/slopes[i] - 1/slopes[i-1])
            
        elif slopes[i] == 0 and slopes[i-1] != 0:
            rhs[i] = c * (-1/hs[i]**2 * dxs[i] - 1/hs[i] * 1/slopes[i-1])
            
        elif slopes[i] != 0 and slopes[i-1] == 0:
            rhs[i] = c * (1/hs[i] * 1/slopes[i] -1/hs[i]**2 * dxs[i-1])
            
    if slopes[N-1] != 0:   
        rhs[N] = c * (1/hs[N]) * (1/slopes[N-1]) - pN
    else:
        rhs[N] = c * (-1/hs[N]**2) * dxs[N-1] - pN
        
    return rhs
        
def make_ps_two(domain, height, vs, N):
    
    ps = np.zeros(domain.Nx)
    slopes = height.slopes
    dxs = height.dxs
    hs = height.hs
    x_peaks = height.x_peaks
    xs = domain.xs
    
    k = 0
    
    cq = -1/2 * vs[-1]
    cu = -6 * domain.eta * domain.U
    
    for i in range(domain.Nx):
        
        if xs[i] > x_peaks[k+1] and k+1 < N:
            k += 1
        
        if slopes[k] != 0:
            dxdh = 1/slopes[k]
            hr = 1/hs[i]
            ps[i] = cq * dxdh * hr**2 + cu * dxdh * hr + vs[k]
        else:
            dx = dxs[k]
            hr = 1/hs[i]
            ps[i] = cq * dx * hr**3 + cu * (-1 * hr**2) * dx + vs[k]
    return ps     
    
        
        
        
        
        
        
        
        
        
        
        
