# -*- coding: utf-8 -*-
"""
Created on Thu May  9 11:59:05 2024

@author: sarah
"""

import numpy as np
    
def linearHeight_solve(height, p0, pN):
    
    h = height.h0
    
    etaU = 6*height.U*height.visc
    
    cq = h**3 * (pN - p0)/(height.xf - height.x0) - etaU*h
    
    ps = np.zeros(height.Nx)
    
    for i in range(height.Nx):
        dx = height.xs[i] - height.x0
        
        ps[i] = (cq * dx / h**3) + (etaU * dx / h**2) + p0
    return ps

def sinsusoidalHeight_solve(height, p0, pN):
    #TODO: ignores boundary pressures
    ps = np.zeros(height.Nx)
    etaU = 6*height.U*height.visc
    
    for i in range(height.Nx):
        h = height.hs[i]
        hx = height.hxs[i]
        ps[i] = -etaU * (h + height.h_mid)/((height.k*height.h_mid)**2*(2 + height.r**2)) * hx / h**2

def stepHeight_solve(height, p0, pN):
    
    ps = np.zeros(height.Nx)

    h_in = height.hs[0]
    h_out = height.hs[-1]

    x0 = height.x0
    xm = height.x_step
    xf = height.xf
    xs = height.xs
    etaU = 6*height.U*height.visc


    m_in_numer = (h_out/h_in)**3 * (p0 - pN)/(xm - xf) - etaU * (h_out - h_in)/h_in**3
    m_in_denom = 1 - (h_out/h_in)**3 * (xm - x0)/(xm - xf)
    
    m_in = m_in_numer/m_in_denom

    m_out = ((xm - x0)*m_in + (p0 - pN))/(xm - xf)

    for i in range(height.Nx):
        if height.xs[i] <= height.x_step:
            ps[i] = m_in * (xs[i] - x0) + p0
        else:
            ps[i] = m_out * (xs[i] - xf) + pN