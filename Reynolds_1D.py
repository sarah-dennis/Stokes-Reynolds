#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 21 09:43:24 2022

@author: sarahdennis
"""

import numpy as np

from matplotlib import pyplot as pp
import _graphics as graph

import domain as dfd

import examples_1D as eg
#------------------------------------------------------------------------------
# Domain 
#------------------------------------------------------------------------------
x0 = 0
xf = 1

BC = 1 # 1:= fixed, 0:= periodic

Nx = 200

U = 1 #lower surface velocity
eta = 1 #viscosity

#------------------------------------------------------------------------------
domain = dfd.Domain(x0, xf, eta, U, Nx, BC)
#------------------------------------------------------------------------------
#height, pressure = eg. corrugated(domain)
#height, pressure = eg.wedge(domain)
#height, pressure = eg.step(domain)
#height, pressure = eg.twoStep(domain)
height, pressure = eg.squareWave(domain)
#------------------------------------------------------------------------------
def reynolds_rhs(domain, height): 
    fs = np.zeros(domain.Nx) #f[x]
    for i in range(domain.Nx):
        fs[i] = 6 * domain.eta * domain.U * height.hxs[i]
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
def solve(domain, height, pressure, RHS=0, FIG=1):

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
    if BC == 0:
        
        # -- set top right corner to D_lower with j = 0
        D[0, domain.Nx-1] = D_lower[0]
        
        # -- set bottom left corner to D_upper at j = N-1 
        D[domain.Nx-1, 0] = D_upper[domain.Nx-1]
        
        # closure: assume sum p'' = 0 
        D[domain.Nx-1, : ] = 1
        fs[domain.Nx-1] = 0
    
    # adjust for fixed pressure boundary ...
    elif BC == 1:
        
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
        
        if FIG == 1: 
            graph.plot_p_h(pressure.ps, ps_numsol, height.hs, domain.xs, height.h_str)
             
            pp.figure()
            title = "Error: Analytic - Numerical Solutions | $N_x=%d$, $dx = %.2f$ \n $%s$"%(domain.Nx, domain.dx, height.h_str)
            graph.plot_2D(pressure.ps-ps_numsol, domain.xs, title, "error")
        
        return inf_norm_err
    
    else: 
        print("Solved Nx=%d"%domain.Nx) 
        
        if FIG==1:
            title = "Numerical Pressure and Height| $N_x=%d$, $dx = %.2f$ \n $%s$ "%(domain.Nx, domain.dx, height.h_str)
            labels = ["$p(x)$", "$h(x)$"]
            graph.plot_2D_multi([ps_numsol, height.hs], domain.xs, title, labels)
        
      
        return ps_numsol

#------------------------------------------------------------------------------
# Convergence
#------------------------------------------------------------------------------   
#BCs = [0: periodic, 1: fixed]

#RHS = [0: reynolds, 1: manf]

def conveg(trials=10, N0=10, BC=1, RHS=0):
    Nx_k = N0
    infNorms = np.zeros(trials)
    dxs = np.zeros(trials)
    dxs_sqr = np.zeros(trials)
    fig=0
    for k in range(trials):
        if k == trials-1: fig=1
        
        domain_k = dfd.Domain(x0, xf, eta, U, Nx_k, BC)
        
        #---- EXAMPLES n----------------------
        #height_k, pressure_k = eg.wedge(domain_k)
        #height_k, pressure_k = eg.corrugated(domain_k)
        height_k, pressure_k = eg.step(domain_k)
        #height_k, pressure_k = eg.twoStep(domain_k)
        #-------------------------------------------------------

        
        infNorms[k] = solve(domain_k, height_k, pressure_k, RHS, 1, FIG=fig)
        
        dxs[k] = domain_k.dx
        dxs_sqr[k] = dxs[k]**2
        
        Nx_k *= 2
    
    pp.figure(trials+1)

    pp.loglog(dxs, dxs, color='b', label='$\mathcal{O}(dx)$')
    pp.loglog(dxs, dxs_sqr, color='g', label='$\mathcal{O}(dx^2)$')
    pp.loglog(dxs, infNorms, color='r', label='$L_{\infty}$ Norm Error')


    
    pp.xlabel('Grid spacing $dx$')

    pp.legend(loc='lower right')
    pp.title("Convergence for %s"%(height_k.h_str), fontweight ="bold")
    