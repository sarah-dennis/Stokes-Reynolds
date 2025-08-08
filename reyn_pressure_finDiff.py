#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 21 09:43:24 2022

@author: sarahdennis
"""
import numpy as np

# Reynolds rhs
def make_rhs(height, U, dP):
    N = height.Nx

    rhs = np.zeros(N)    

    rhs = 6 * U * height.hxs #*visc
    
    
    p0=-dP
    pN = 0
    rhs[0] = p0
    rhs[N-1] = pN
    
    return rhs

def make_mat(height):
    N = height.Nx    
    D_lower = np.zeros(N)
    D_center = np.zeros(N)
    D_upper = np.zeros(N)
    
    for i in range(N): #space: xs[i]
    
        # Find h(t,x) at x = [xs[i-1], xs[i], xs[i+1]] = [hl, hc, hr]
        hl = height.hs[(i-1) % N]
        hc = height.hs[i % N]   
        hr = height.hs[(i+1) % N]

        #P(i) central diagonal

        D_center[i] = -(hr**3 + 2*hc**3 + hl**3) 
        
        #P(i+1) upper diagonal
        D_upper[i] = (hr**3 + hc**3)
        
        #P(i-1) lower diagonal
        D_lower[i] = (hl**3 + hc**3)

        
    # combine as upper, middle, lower diagonals
    D = np.diagflat(D_center) + np.diagflat(D_lower[1:N], -1) + np.diagflat(D_upper[0:N-1], 1)
    D /= (2*height.dx**2)
    
    # -- set top row D to [1, 0, ...] and rhs[0] = p(x0)
    D[0,0] = 1
    D[0,1] = 0
    # -- set bottom row D to [ ... , 0, 1] and rhs[Nx-1] = p(xf)
    D[N-1,N-1] = 1
    D[N-1,N-2] = 0

    return D
        
def make_rhs_Q(height, U, Q):
    N = height.Nx

    rhs = np.zeros(N)    

    rhs = 6 * U * height.hxs  #*visc
    h0 = height.hs[0]
    
    dp0 = -12*Q/h0**3 + 6*U/h0**2
    
    rhs[0] = dp0#*visc
    rhs[N-1] = 0
    
    return rhs


def make_mat_Q(height):
    N = height.Nx    
    D_lower = np.zeros(N)
    D_center = np.zeros(N)
    D_upper = np.zeros(N)
    
    for i in range(N): #space: xs[i]
    
        # Find h(t,x) at x = [xs[i-1], xs[i], xs[i+1]] = [hl, hc, hr]
        hl = height.hs[(i-1) % N]
        hc = height.hs[i % N]   
        hr = height.hs[(i+1) % N]

        #P(i) central diagonal

        D_center[i] = -(hr**3 + 2*hc**3 + hl**3) 
        
        #P(i+1) upper diagonal
        D_upper[i] = (hr**3 + hc**3)
        
        #P(i-1) lower diagonal
        D_lower[i] = (hl**3 + hc**3)

        
    # combine as upper, middle, lower diagonals
    D = np.diagflat(D_center) + np.diagflat(D_lower[1:N], -1) + np.diagflat(D_upper[0:N-1], 1)
    D /= (2*height.dx**2)
    
    # -- set top row D to [-3/2, 2, -1/2 ...] and rhs[0] = p(x0)
    D[0,0] = -3/(2*height.dx)
    D[0,1] = 4/(2*height.dx)
    D[0,2] = -1/(2*height.dx)
    # -- set bottom row D to [ ... , 0, 1] and rhs[Nx-1] = p(xf)
    D[N-1,N-1] = 1
    D[N-1,N-2] = 0

    return D


