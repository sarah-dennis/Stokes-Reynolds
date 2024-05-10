# -*- coding: utf-8 -*-
"""
Created on Thu May  9 11:59:05 2024

@author: sarah
"""

import numpy as np
from pressures import P_Solver

class Solver_Sinusoidal(P_Solver):
    def __init__(self, height, p0, pN):
        super().__init__(height, p0, pN, self.solve)

    def solve(self):
        #TODO: ignores boundary pressures
        # p0 = self.p0 = 0
        # pN = self.pN = 0
        Nx = self.height.Nx
        hs = self.height.hs
        hxs = self.height.hxs
        h_mid = self.height.h_mid
        k = self.height.k
        r = self.height.r
        etaU = 6*self.height.U*self.height.visc
        
        ps = np.zeros(Nx)
        for i in range(Nx):
            h = hs[i]
            hx = hxs[i]
            ps[i] = -etaU * (h + h_mid) * h**2 / (hx * (k*h_mid)**2 * (2 + r**2)) 
        return ps
    
class Solver_Constant(P_Solver):
    def __init__(self, height, p0, pN):
        super().__init__(height, p0, pN, self.solve)
    
    def solve(self):
        p0 = self.p0
        pN = self.pN
        x0 = self.height.x0
        xf = self.height.xf
        xs = self.height.xs
        h0 = self.height.h0
        Nx = self.height.Nx
        
        etaU = 6*self.height.U*self.height.visc
        
        cq = h0**3 * (pN - p0)/(xf - x0) - etaU*h0
        
        ps = np.zeros(Nx)

        for i in range(Nx):
            dx = xs[i] - x0
            ps[i] = (cq * dx / h0**3) + (etaU * dx / h0**2) + p0
        return ps

    
class Solver_Step(P_Solver):
    def __init__(self, height, p0, pN):
        super().__init__(height, p0, pN, self.solve) 
    
    def solve(self):
        p0 = self.p0
        pN = self.pN
        Nx = self.height.Nx
        
        ps = np.zeros(Nx)
    
        h_in = self.height.hs[0]
        h_out = self.height.hs[-1]
    
        x0 = self.height.x0
        xm = self.height.x_step
        xf = self.height.xf
        xs = self.height.xs
        etaU = 6*self.height.U*self.height.visc
    
    
        m_in_numer = (h_out/h_in)**3 * (p0 - pN)/(xm - xf) - etaU * (h_out - h_in)/h_in**3
        m_in_denom = 1 - (h_out/h_in)**3 * (xm - x0)/(xm - xf)
        
        m_in = m_in_numer/m_in_denom
    
        m_out = ((xm - x0)*m_in + (p0 - pN))/(xm - xf)
    
        for i in range(Nx):
            if xs[i] <= xm:
                ps[i] = m_in * (xs[i] - x0) + p0
            else:
                ps[i] = m_out * (xs[i] - xf) + pN
        return ps