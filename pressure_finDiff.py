#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 21 09:43:24 2022

@author: sarahdennis
"""
import numpy as np
import time

#------------------------------------------------------------------------------
# Numerical solution
#------------------------------------------------------------------------------
# RHS = [0: reynolds, 1: exact]
# FIG = [0: two plots (Numerical and analytic pressure, and Error), 1: one plot only finDiff pressure]

def solve(height, p0, pN):
    N = height.Nx
    
    t0 = time.time()
    
    # initialise RHS = 6 eta U hx
    fs = np.zeros(N)
    for i in range(N):
        fs[i] = 6 * height.eta * height.U * height.hxs[i]

    # initilise diagonals of differnce matrix
    D_lower = np.ones(N)
    D_center = np.ones(N)
    D_upper = np.ones(N)
    
    for i in range(N): #space: xs[i]
    
        # Find h(t,x) at x = [xs[i-1], xs[i], xs[i+1]] = [hl, hc, hr]
        hl = height.hs[(i-1) % N]
        hc = height.hs[i % N]   
        hr = height.hs[(i+1) % N]

        #P(i) central diagonal
        D_center[i] = -(hr**3 + 2*hc**3 + hl**3)/(2*(height.dx**2)) 
        
        #P(i+1) upper diagonal
        D_upper[i] = (hr**3 + hc**3)/(2*(height.dx**2))
        
        #P(i-1) lower diagonal
        D_lower[i] = (hl**3 + hc**3)/(2*(height.dx**2))
        
        
    # combine as upper, middle, lower diagonals
    D = np.diagflat(D_center) + np.diagflat(D_lower[1:N], -1) + np.diagflat(D_upper[0:N-1], 1)
        
    # -- set top row D to [1, 0, ...] and f[0] = p_inlet
    D[0,0] = 1
    D[0,1] = 0
    fs[0] = p0
        
    # -- set bottom row D to [ ... , 0, 1] and f[Nx-1] = p_outlet
    D[N-1,N-1] = 1
    D[N-1,N-2] = 0
    fs[N-1] = pN
        
    # solve for p
    ps_numsol = np.linalg.solve(D, fs)
    tf = time.time()
    
    print("Solved Numerically for Nx=%d in %.5f \n"%(N, tf-t0))
    return ps_numsol





