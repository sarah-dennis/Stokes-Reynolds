#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 29 13:14:58 2024

@author: sarahdennis
"""

        # Xs = height.x_peaks #[x0, x1, ..., xN]
        # Hs = height.h_peaks #[h0, h1, ..., hN]
        # slopes = height.slopes #[(i, i+1): i = 0, ... N-1]
        
        # N = len(Xs)
        
        # rhs = st.make_rhs(N, Xs, Hs, slopes, p0, pN)
        
        # st_linOp = st.stLinOp(N, Xs, Hs, slopes)
        
        # Cs = spsolve(st_linOp, rhs) #N+1 : [Cp0, Cp1, ..., Cp_N-1, Cq]
        
        # ps = st.make_ps(domain, height, Cs)
        
from scipy.sparse.linalg import LinearOperator 
import numpy as np

class stLinOp(LinearOperator):
    def __init__(self, N, xs, hs, slopes):
        self.shape = (N,N)
        self.dtype = np.dtype('f8')
        self.N = N
        self.xs = xs
        self.hs = hs
        self.hs_sqrd = np.array([h**2 for h in hs])
        
        self.slopes = slopes
        
        
        self.rcp_slopes = self.make_slope_reciprocals(N, xs, hs, slopes)
        
        self.mv = np.zeros(N)
        
    def make_slope_reciprocals(N, xs, hs, slopes):
        
        

    def _matvec(self, c):
        xs = self.xs
        hs = self.hs
        hhs = self.hs_sqrd
        N = self.N
        
        mc = np.zeros(N)
        
        mc[0] = c[0] - 0.5
        
        for i in range(N):
            