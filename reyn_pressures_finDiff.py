#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 21 09:43:24 2022

@author: sarahdennis
"""
import numpy as np


def make_rhs(height):
    N = height.Nx
    c = 6 * height.visc * height.U* height.dx

    rhs = np.zeros(N)    
    rhs[0] = height.p0
    for i in range(1,N-1):
        rhs[i] = c*(height.hs[i+1]-height.hs[i-1])
    rhs[N-1] = height.pN
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
        D_center[i] = -(hr**3 + 2*hc**3 + hl**3) #/(2*(height.dx**2)) 
        
        #P(i+1) upper diagonal
        D_upper[i] = (hr**3 + hc**3) #/(2*(height.dx**2))
        
        #P(i-1) lower diagonal
        D_lower[i] = (hl**3 + hc**3) #/(2*(height.dx**2))
        
        
    # combine as upper, middle, lower diagonals
    D = np.diagflat(D_center) + np.diagflat(D_lower[1:N], -1) + np.diagflat(D_upper[0:N-1], 1)
        
    # -- set top row D to [1, 0, ...] and f[0] = p_inlet

    D[0,0] = 1
    D[0,1] = 0
    # -- set bottom row D to [ ... , 0, 1] and f[Nx-1] = p_outlet
    D[N-1,N-1] = 1
    D[N-1,N-2] = 0

    return D
        





