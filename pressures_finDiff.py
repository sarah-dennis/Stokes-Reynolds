#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 21 09:43:24 2022

@author: sarahdennis
"""
import numpy as np

from pressures import P_Solver
#------------------------------------------------------------------------------
# Numerical solution
#------------------------------------------------------------------------------

class Solver_finDiff(P_Solver):
    def __init__(self, height, p0, pN):
        p_str = "Reynolds Finite Difference ($N_x = %d$)"%height.Nx
        
        self.rhs = make_rhs(height, p0, pN)
        self.fdMat = make_fdMat(height, p0, pN)
        
        super().__init__(height, p0, pN, self.solve, p_str)
                
    def solve(self):
        ps = np.linalg.solve(self.fdMat, self.rhs)
        return ps
    
def make_rhs(height, p0, pN):
    N = height.Nx
    c = 6 * height.visc * height.U

    rhs = np.zeros(N)
    for i in range(N):
        rhs[i] = c * height.hxs[i]
        
    rhs[0] = p0
    rhs[N-1] = pN
    return rhs

def make_fdMat(height, p0, pN):
    N = height.Nx    

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
        
    # -- set bottom row D to [ ... , 0, 1] and f[Nx-1] = p_outlet
    D[N-1,N-1] = 1
    D[N-1,N-2] = 0
    
    return D
        





