# -*- coding: utf-8 -*-
"""
Created on Tue Jul 30 15:17:44 2024

@author: sarah
"""

import numpy as np
import time
        
from scipy.sparse.linalg import LinearOperator 
from scipy.sparse.linalg import gmres
import reyn_boundary as bc


def gmres_solve(height, BC):
    t0 = time.time()
    rhs = make_rhs(height, BC)
    linOp = pwlLinOp(height,BC)
    sol_coefs, exit_code = gmres(linOp, rhs, tol=1e-10)
        
    if exit_code != 0:
        raise Exception('gmres did not converge')
   
    tf = time.time()
    print('gmres time: ', tf-t0)
    ps_1D = make_ps(height, BC, sol_coefs)
    tF =time.time()
    # print('gmres total time: ', tF-t0)
    return ps_1D

def make_ps(height, BC, cs):
    t0=time.time()
    slopes = height.slopes
    hs = height.hs
    x_peaks = height.x_peaks
    xs = height.xs
    
    ps = np.zeros(height.Nx)
    k = 0
    cq = cs[-1]
    #flux = -cq/12 #/visc
    cu = 6 * BC.U #*visc
   
    for i in range(height.Nx):
        
        if xs[i] > x_peaks[k+1]:
            k += 1
        
        if slopes[k] != 0:
            dhdx = slopes[k]
            h = hs[i]
            ps[i] = -(cq/2 *h**-2 + cu/h)/dhdx + cs[k]
        else: 
            dx = xs[i] - x_peaks[k]
            h = hs[i]
            ps[i] = dx* (cq * h**-3 + cu * h**-2) + cs[k]
    tf=time.time()
    # print('make ps time:', tf-t0)
    return ps


def make_rhs(height, BC):
    N = height.N_regions
    hs = height.h_peaks
    slopes = height.slopes
    widths = height.widths
    
    rhs = np.zeros(N+1)
    c = 6*BC.U #*visc
    if isinstance(BC, bc.Fixed):
        if slopes[0] != 0:
            rhs[0] = c / (hs[0,1] * slopes[0]) + BC.p0
        else: 
            rhs[0] = BC.p0 
        
                    
        if slopes[N-1] != 0:
            rhs[N] = c / (hs[N,0]*slopes[N-1]) + BC.pN
        else:
            rhs[N] = -c * hs[N,0]**-2 * widths[N-1] + BC.pN
            
    elif isinstance(BC, bc.Mixed):
        rhs[0] = -12 * BC.Q
        if slopes[N-1] != 0:
            rhs[N] = c / (hs[N,0]*slopes[N-1]) + BC.pN
        else:
            rhs[N] = -c * hs[N,0]**-2 * widths[N-1] + BC.pN
    
    for i in range(1, N-1):
        
        if slopes[i] != 0 and slopes[i-1] != 0:
            rhs[i] = c *(1/(hs[i,1] * slopes[i]) - 1/(hs[i,0] * slopes[i-1]))
            
        elif slopes[i] != 0 and slopes[i-1] == 0:
            rhs[i] = c * (1/(hs[i,1] * slopes[i]) + hs[i,0]**-2 * widths[i-1])

        elif slopes[i] == 0 and slopes[i-1] != 0:
            rhs[i] = -c * ( 1/(hs[i,0] * slopes[i-1]))

        else: 
            rhs[i] = c * hs[i,0]**-2 * widths[i-1]
  
    return rhs

class pwlLinOp(LinearOperator):
    def __init__(self, height,BC):

        self.N = height.N_regions
        self.h_peaks = height.h_peaks
        self.widths = height.widths
        self.slopes = height.slopes
        self.BC = BC
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
        
        if isinstance(self.BC, bc.Fixed):
            if slopes[0] != 0:
                mv[0] = v[0] - cq/(2 * hs[0,1]**2 * slopes[0])
            else:
                mv[0] = v[0]
                  
            if slopes[N-1] != 0:
                mv[N] = v[N-1] - cq/(2 * hs[N,0]**2 * slopes[N-1])
            else:
                mv[N] = v[N-1] + cq * hs[N,0]**-3 * widths[N-1]
        
        elif isinstance(self.BC, bc.Mixed):
            mv[0] = cq
            if slopes[N-1] != 0:
                mv[N] = v[N-1] - cq/(2 * hs[N,0]**2 * slopes[N-1])
            else:
                mv[N] = v[N-1] + cq * hs[N,0]**-3 * widths[N-1]
 
        for i in range(1, N):
            if slopes[i] != 0 and slopes[i-1] != 0:
                mv[i]= -v[i-1] + v[i] - cq/2 * (1/(hs[i,1]**2 * slopes[i]) - 1/(hs[i,0]**2 * slopes[i-1]))
            
            elif slopes[i] != 0 and slopes[i-1] == 0:
                mv[i] = -v[i-1] + v[i] - cq * ( 1/(2* hs[i,1]**2 * slopes[i]) + hs[i,0]**-3 * widths[i-1])

            elif slopes[i] == 0 and slopes[i-1] != 0:
                mv[i] = -v[i-1] + v[i] + cq * (1/(2*hs[i,0]**2 * slopes[i-1]))
                
            else:
                mv[i] = -v[i-1] + v[i] - cq * hs[i,0]**-3 * widths[i-1]
           
        return mv

