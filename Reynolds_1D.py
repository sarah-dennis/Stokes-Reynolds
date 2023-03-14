#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 21 09:43:24 2022

@author: sarahdennis
"""
import numpy as np

import _graphics as graph

import domain as dfd
import examples_1D as eg

#------------------------------------------------------------------------------
# Domain & Discretization
#------------------------------------------------------------------------------
x0 = 0      # left boundary 
xf = 1      # right boundary

Nx = 8000   # Number of grid points

BC = "fixed"      
#BC = "periodic"

U = 1      # lower surface velocity
eta = 1     # flouid viscosity

domain = dfd.Domain(x0, xf, eta, U, Nx, BC)


#------------------------------------------------------------------------------
# Height & Pressure
#------------------------------------------------------------------------------
p0 = 0
pN = 0

#height, pressure = eg.corrugated(domain, p0, pN)
#height, pressure = eg.wedge(domain, p0, pN)
#height, pressursole = eg.step(domain, p0, pN)
#height, pressure = eg.twoStep(domain, p0, pN)
#height, pressure = eg.squareWave(domain, p0, pN)
#height, pressure = eg.dimple(domain, p0, pN)
height, pressure = eg.flat(domain, p0, pN)


#------------------------------------------------------------------------------
# RHS
#------------------------------------------------------------------------------
def reynolds_rhs(domain, height): 
    fs = np.zeros(domain.Nx) #f[x]
    for i in range(domain.Nx):
        fs[i] = 6 * domain.eta * domain.U * height.hxs[i]
        #TODO ^^ express using Re
    return fs

def manf_rhs(domain, height, pressure):
    fs = np.zeros(domain.Nx) #f[x]
    for i in range(domain.Nx):
        fs[i] = (height.hs[i] ** 3) * pressure.pxxs[i] + 3 * height.hs[i]**2 * height.hxs[i] * pressure.pxs[i]
    return fs

#------------------------------------------------------------------------------
# Numerical solution
#------------------------------------------------------------------------------
# RHS = [0: reynolds, 1: exact]
def solve(domain=domain, height=height, pressure=pressure, RHS=0, FIG=1):

    # Reynolds RHS = 6 eta U hx
    if RHS == 0: 
        fs = reynolds_rhs(domain, height)
    
    # ManfSol RHS
    elif RHS == 1: 
         fs = manf_rhs(domain, height, pressure)

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
    if BC == "periodic":
        
        # -- set top right corner to D_lower with j = 0
        D[0, domain.Nx-1] = D_lower[0]
        
        # -- set bottom left corner to D_upper at j = N-1 
        D[domain.Nx-1, 0] = D_upper[domain.Nx-1]
        
        # closure: assume sum p'' = 0 
        D[domain.Nx-1, : ] = 1
        fs[domain.Nx-1] = 0
    
    # adjust for fixed pressure boundary ...
    elif BC == "fixed":
        
        # -- set top row D to [1, 0, ...] and f[0] = p_inlet
        D[0,0] = 1
        D[0,1] = 0
        fs[0] = pressure.p0
        
        # -- set bottom row D to [ ... , 0, 1] and f[Nx-1] = p_outlet
        D[domain.Nx-1,domain.Nx-1] = 1
        D[domain.Nx-1,domain.Nx-2] = 0
        fs[domain.Nx-1] = pressure.pf

    # solve for p
    ps_numsol = np.linalg.solve(D, fs)
    
    # Plotting and error 
    inf_norm_err = np.max(np.abs(pressure.ps-ps_numsol))
    print ("Solved Nx=%d with error %0.5f"%(domain.Nx, inf_norm_err))
    #  unknown exact pressure: ps = [0, ..., 0]
        
    if FIG == 1: 
        p_h_title = "Numerical Reynolds and Exact Pressure for %s"%height.h_str
        graph.plot_pN_pE_h(pressure.ps, ps_numsol, height.hs, domain.xs, p_h_title)

        err_title = "Error: %s | $N_x=%d$, $dx = %.3f$"%(height.h_str, domain.Nx, domain.dx)
        graph.plot_2D(pressure.ps-ps_numsol, domain.xs, err_title, "error", "dx")
    elif FIG == 2:
        p_h_title = "Numerical Reynolds for %s"%height.h_str
        graph.plot_p_h(ps_numsol, height.hs, domain.xs, p_h_title)

    return inf_norm_err, ps_numsol