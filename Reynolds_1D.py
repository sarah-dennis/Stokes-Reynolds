#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 21 09:43:24 2022

@author: sarahdennis
"""

import numpy as np
from matplotlib import pyplot as pp

#------------------------------------------------------------------------------
# Domain & Height
#------------------------------------------------------------------------------
# Grid sizing Nx given at run time

# width
xa, xb = (0, 2*np.pi)

# Height
# y = h(x)
# -- option 1: Sinusoidal height 
# h0, delta, k = (0.5, 0.1, 4) # delta < h0 so height positive
# str_h = "%0.1f + %0.1f \cos(%d x)"%(h0, delta, k) #for graph title

# def h(x):
#     return h0 + delta * np.sin(k*x)

# def h_dx(x):
#     return delta * k * np.cos(k*x)

# -- option 2: Wedge height 
h0, hf = (0.01, 0.001) #h0 = inital height, hf: final height
delta = (hf - h0)/(xb - xa) # appropriate slope

def h(x):
    return delta * (x - xa) + h0

def h_dx(x):
    return delta
str_h = "%0.1f + %0.2f (x - %d)"%(h0, delta, xa) #for graph title
                
#------------------------------------------------------------------------------
# Manf. solution & RHS
#------------------------------------------------------------------------------
# (h^3 u')' = f(h)

alpha = 2
def p(x):
    return -np.sin(alpha*x)

str_p = "-\sin(%dx)"%(alpha) #for graph title

def p_dx(x):
    return -alpha * np.cos(alpha*x)

def p_dxx(x): 
    return alpha**2 * np.sin(alpha*x)

def f(x): # = h^3 u'' + (h^3)' u' 
    return (h(x) ** 3) * p_dxx(x) + 3 * h(x)**2 * h_dx(x) * p_dx(x)


#------------------------------------------------------------------------------
# numerical solution
#------------------------------------------------------------------------------   
def solve(Nx=100, figrs=1, BC=0):
    
    dx = (xb - xa)/(Nx-1)

    xs = [xa + i*dx for i in range(Nx)]
    
    # construct height on grid
    hs = np.zeros(Nx) #h[t][x]
    for i in range(Nx):
            hs[i] = h(xs[i])
    
    # construct f on grid
    f_n = np.zeros(Nx) #f[x]
    for i in range(Nx):
        f_n[i] = f(xs[i])
        
    # construct exact solution on grid
    p_exact = np.zeros(Nx)
    for i in range(Nx):
        p_exact[i] = p(xs[i])
                
    inf_norm_err = 0
    
    # initilise diagonals of differnce matrix
    D_lower = np.ones(Nx)
    D_center = np.ones(Nx)
    D_upper = np.ones(Nx)
    
    for i in range(Nx): #space: xs[i]
    
        # Find h(t,x) at x = [xs[i-1], xs[i], xs[i+1]] = [hl, hc, hr]
        hl = hs[(i-1) % Nx] 
        hc = hs[i % Nx]       
        hr = hs[(i+1) % Nx] 
        
        # option 1: h^3 = [(h + h)/2]^3
        #P(i) central diagonal
        # D_center[i] = -(hr**3 + 3*hc * (hr**2 + hr*hc + hc*hl + hl**2) + 2*hc**3 + hl**3)/(8*(dx**2))
        
        # #P(i+1) upper diagonal
        # D_upper[i] = (hr**3 + 3*hc*hr * (hc + hr) + hc**3)/(8*(dx**2))
    
        # #P(i-1) lower diagonal
        # D_lower[i] = (hl**3 + 3*hc*hl * (hc + hl) + hc**3)/(8*(dx**2))
        
        #option 2: h^3 = (h^3 + h^3)/2
        #P(i) central diagonal
        D_center[i] = -(hr**3 + 2*hc**3 + hl**3)/(2*(dx**2)) 
        
        #P(i+1) upper diagonal
        D_upper[i] = (hr**3 + hc**3)/(2*(dx**2))
        
        #P(i-1) lower diagonal
        D_lower[i] = (hl**3 + + hc**3)/(2*(dx**2))
        
        
    # combine as upper, middle, lower diagonals
    D = np.diagflat(D_center) + np.diagflat(D_lower[1:Nx], -1) + np.diagflat(D_upper[0:Nx-1], 1)

    
    # adjust for periodic boundary...
    if BC == 0:
        
        # -- set top right corner to D_lower with j = 0
        D[0, Nx-1] = D_lower[0]
        
        # -- set bottom left corner to D_upper at j = N-1 
        D[Nx-1, 0] = D_upper[Nx-1]
        
        # closure: assume sum p'' = 0 
        D[Nx-1, : ] = 1
        f_n[Nx-1] = 0
    
    # adjust for fixed pressure boundary
    elif BC == 1:
        
        # -- set top row D to [1, 0, ...] and f[0] = p_inlet
        D[0,0] = 1
        D[0,1] = 0
        f_n[0] = p(xs[0])
        
        # -- set bottom row D to [ ... , 0, 1] and f[Nx-1] = p_outlet
        D[Nx-1,Nx-1] = 1
        D[Nx-1,Nx-2] = 0
        f_n[Nx-1] = p(xs[Nx-1])

    # solve for p
    p_n = np.linalg.solve(D, f_n)
    #return p_n
    
    inf_norm_err = np.max(np.abs(np.subtract(p_exact,p_n)))
    print ("Solved Nx=%d with error %0.5f"%(Nx, inf_norm_err))    

    if figrs: 
        
        pp.figure()
        
        pp.plot(xs, p_exact, label="$p$", color='r')
        
        pp.plot(xs, p_n, label="$p_N$", color='b');
        
        #plot h^3 u'
        # hs_cube = [h**3 for h in hs
        # p_dxs = [p_dx(x) for x in xs]
        # inner = [hs_cube[i] * p_dxs[i] for i in range(Nx)]
        # pp.plot(xs, inner, label="$h^3 p'$", color='black');
        
        #plot f
        # pp.plot(xs, f_n, label="f", color='g')
        
        #plot h
        pp.plot(xs, hs, label="h", color='g')
        
        pp.xlabel('x')
        pp.ylabel('p(x)')
        pp.legend()
        
        pp.title("$(h^3p')'= f$ | $p = %s$ | $h=%s$ | $N_x=%d$"%(str_p, str_h, Nx))

    return inf_norm_err

#------------------------------------------------------------------------------
# Convergence
#------------------------------------------------------------------------------   
def conveg(trials=10, N0=5, BC=1, figrs=0):
    N = N0
    infNorms = np.zeros(trials)
    dxs = np.zeros(trials)
    dxs_sqr = np.zeros(trials)
    
    for i in range(trials):
        
        infNorms[i] = solve(N, figrs, BC)
        dxs[i] = ((xb - xa) / N)
        dxs_sqr[i] = dxs[i]**2
        
        N *= 2
    
    pp.figure(trials+1)
    pp.loglog(dxs, dxs_sqr, color='r', label='$dx^2$')
    pp.loglog(dxs, infNorms, color='b', label='Linf Error')
    
    pp.xlabel('dx')

    pp.legend()
    pp.title('$N_x$=%i to $N_x$=%i with %i trials'%(N0, N/2, trials))



