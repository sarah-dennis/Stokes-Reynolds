#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 21 09:43:24 2022

@author: sarahdennis
"""

import numpy as np
from matplotlib import pyplot as pp

#------------------------------------------------------------------------------
# Domain
#------------------------------------------------------------------------------
xa, xb = (0, 2*np.pi)

#------------------------------------------------------------------------------
# manuf. solution
#------------------------------------------------------------------------------
alpha = 4
def p(x):
    return -np.sin(alpha*x)
    
str_p = "-\sin(%dx)"%(alpha) #for graph title

def f(x):  #u'' = f
    return alpha**2 * np.sin(alpha*x)
 
    #return alpha**2 * np.cos(alpha*x)     

#------------------------------------------------------------------------------
# Discretization  
#------------------------------------------------------------------------------   
def solve(Nx=100, fig_k=0):
    print(f'solving  N = {Nx}')
    dx = (xb - xa)/(Nx)

    xs = [xa + i*dx for i in range(Nx)]
    
    #construct f on grid
    f_n = np.zeros((Nx, 1))
    for i in range (Nx):
        f_n[i]= f(xs[i])
        
    #construct exact solution on grid
    p_exact = np.zeros((Nx, 1))
    for i in range(Nx):
        p_exact[i] = p(xs[i])
    
    # construct diag(1, -2, 1)
    D_mid = np.diagflat(np.ones(Nx)*-2)
    D_upper = np.diagflat(np.ones(Nx-1)*1, 1)
    D_lower = np.diagflat(np.ones(Nx-1)*1, -1)
    D =  D_mid + D_lower + D_upper

    # adjust for periodic boundary
    D[0,Nx-1] = 1
    D[Nx-1, 0] = 1
    
    D = D / (dx**2)
    
    #assume sum pi = 0
    D[Nx-1, 0:Nx] = 1
    f_n[Nx-1] = 0
    
    #solve for u
    p_n = np.linalg.solve(D, f_n) 
    
    #plotting
    if fig_k != 0:
        pp.figure(fig_k)
        pp.plot(xs,p_exact[0:Nx], label='$p$', color='r')
        pp.plot(xs ,p_n[0:Nx], label='$p_N$', color='b');
        
        pp.xlabel('x')
        pp.ylabel('p(x)')
        pp.legend()
        
        pp.title('Nx=%i, dx=%.4f'%(Nx,dx))
        
        pp.title("$p''= f$ | $p = %s$ | $N_x=%d$"%(str_p, Nx))
        
    
    # error     
    infNorm = np.max(np.abs(p_exact-p_n) );
    print(f'Error {infNorm:.5f}')
    return infNorm

#------------------------------------------------------------------------------
# Convergence
#------------------------------------------------------------------------------   
def conveg(trials=10, Nx0=6):
    Nx = Nx0
    infNorms = np.zeros(trials)
    dxs = np.zeros(trials)
    dxs_p = np.zeros(trials)
    
    for i in range(trials):
        
        infNorms[i] = solve(Nx, i)
        dxs[i] = (xb - xa) / (Nx-1)
        dxs_p[i] = dxs[i] ** 2
        
        Nx *=2
    
    pp.figure(trials)
    pp.loglog(dxs, dxs_p, color='r', label='dx^2')
    pp.loglog(dxs, infNorms, color='b', label='Linf Error')
    
    pp.xlabel('dx')

    pp.legend()
    pp.title('Nx0=%i to Nx=%i with %i trials'%(Nx0,Nx/2,trials))


















