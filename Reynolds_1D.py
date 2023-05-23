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
# FIG = [0: two plots (Numerical and analytic pressure, and Error), 1: one plot only Reynolds pressure]

def solve(domain, height, p0, pN):
    t0 = time.time()
    
    # initialise RHS = 6 eta U hx
    fs = np.zeros(domain.Nx)
    for i in range(domain.Nx):
        fs[i] = 6 * domain.eta * domain.U * height.hxs[i]

    # initilise diagonals of differnce matrix
    D_lower = np.ones(domain.Nx)
    D_center = np.ones(domain.Nx)
    D_upper = np.ones(domain.Nx)
    
    for i in range(domain.Nx): #space: xs[i]
    
        # Find h(t,x) at x = [xs[i-1], xs[i], xs[i+1]] = [hl, hc, hr]
        hl = height.hs[(i-1) % domain.Nx]
        hc = height.hs[i % domain.Nx]   
        hr = height.hs[(i+1) % domain.Nx]

        #P(i) central diagonal
        D_center[i] = -(hr**3 + 2*hc**3 + hl**3)/(2*(domain.dx**2)) 
        
        #P(i+1) upper diagonal
        D_upper[i] = (hr**3 + hc**3)/(2*(domain.dx**2))
        
        #P(i-1) lower diagonal
        D_lower[i] = (hl**3 + hc**3)/(2*(domain.dx**2))
        
        
    # combine as upper, middle, lower diagonals
    D = np.diagflat(D_center) + np.diagflat(D_lower[1:domain.Nx], -1) + np.diagflat(D_upper[0:domain.Nx-1], 1)

    # adjust for periodic boundary...
    if domain.BC == "periodic":
        
        # -- set top right corner to D_lower with j = 0
        D[0, domain.Nx-1] = D_lower[0]
        
        # -- set bottom left corner to D_upper at j = N-1 
        D[domain.Nx-1, 0] = D_upper[domain.Nx-1]
        
        # closure: assume sum p'' = 0 
        D[domain.Nx-1, : ] = 1
        fs[domain.Nx-1] = 0
    
    # adjust for fixed pressure boundary ...
    elif domain.BC == "fixed":
        
        # -- set top row D to [1, 0, ...] and f[0] = p_inlet
        D[0,0] = 1
        D[0,1] = 0
        fs[0] = p0
        
        # -- set bottom row D to [ ... , 0, 1] and f[Nx-1] = p_outlet
        D[domain.Nx-1,domain.Nx-1] = 1
        D[domain.Nx-1,domain.Nx-2] = 0
        fs[domain.Nx-1] = pN
        
    # solve for p
    ps_numsol = np.linalg.solve(D, fs)
    tf = time.time()
    
    print ("Solved Numerically for Nx=%d in %.5f \n"%(domain.Nx, tf-t0))
    return ps_numsol





