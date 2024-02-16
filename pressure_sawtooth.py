#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 15 11:52:42 2024

@author: sarahdennis
"""
import numpy as np
from scipy.sparse.linalg import LinearOperator

class sawtoothLinOp(LinearOperator):
    #xs:= x where stitching occurs
    
    #hs:= height at xi

    def __init__(self, n, xs, hs, U):
        self.n = n

        self.xs = xs #len n+1
        self.hs = hs #len n+1
        self.U = U
        self.visc=1
        self.a, self.b = self.makeABs()
        
        self.shape = (n, n)
        self.dtype = np.dtype('f8')
        self.mv = np.zeros(n)

    
    def makeABs(self):
        a = np.array(self.n)
        b = np.array(self.n)
        
        c = -6*self.U*self.visc
        
        for i in range(self.n):
            xi = self.xs[i]
            xii =self.xs[i+1]
            hi = self.hs[i]
            hii = self.hs[i+1]
            
            a[i] = c * (xii - xi)/(hii - hi) * (1/hii + 1/hi)
        
            b[i] = (1/hii - 1/hi)/2
        return a, b

    def makeRHS(self):
        

# v = [p1, p2, ..., pN-1, hpmax]
    def _matvec(self, v):
        n = self.n
        mv = np.zeros(n)
        
        for i in range(self.n):
            ab = self.a[i]*self.b[i]
            if i == 0:
                mv[i] = -v[0] + ab * v[n-1]
            elif i == n-1:
                mv[i] = v[n-2] + ab * v[n-1] 
            else:
                mv[i]= -v[i] + v[i+1] + ab * v[n-1]
                
        return mv