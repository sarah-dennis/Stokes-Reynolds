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
# width
xa, xb = (0, 4*np.pi)

#time
t0, tf = (0, 2*np.pi)

# height h(t,x) = h(x)  # delta < h0 so height positive
h0, delta, k = (0.15, 0.1, 4)
str_h = "%0.1f + %0.1f \cos(%d x)"%(h0, delta, k) #for graph title

def h(t, x):
    return h0 + delta * np.sin(k*x)

def h_dx(t,x):
    return delta * k * np.cos(k*x)
                
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

def f(t,x): # = h^3 u'' + (h^3)' u' 
    return (h(t,x) ** 3) * p_dxx(x) + 3 * h(t,x)**2 * h_dx(t,x) * p_dx(x)


#------------------------------------------------------------------------------
# numerical solution
#------------------------------------------------------------------------------   
def solve(Nx=100, Nt = 1, figrs=1):
    
    dx = (xb - xa)/Nx
    dt = (tf-t0)/Nt

    xs = [xa + i*dx for i in range(Nx)]
    ts = [t0 + i*dt for i in range(Nt)]
    
    # construct height on grid
    hs = np.zeros((Nt,Nx)) #h[t][x]
    for i in range(Nt): 
        for j in range(Nx):
            hs[i][j] = h(ts[i], xs[j])
    
    # construct f on grid
    f_n = np.zeros((Nt, Nx)) #f[t][x]
    for i in range (Nt):
        for j in range(Nx):
            f_n[i][j] = f(ts[i], xs[j])
        
    # construct exact solution on grid
    p_exact = np.zeros(Nx)
    for i in range(Nx):
        p_exact[i] = p(xs[i])
                
             
    # At each time t = ts[i]
    
    # construct finite difference matrix with h(t)
    inf_norms = np.zeros(Nt)  
    
    for i in range(Nt): #time: ts[i]
    
        # initilise diagonals of differnce matrix
        D_lower = np.ones(Nx-1)
        D_center = np.ones(Nx)
        D_upper = np.ones(Nx-1)
        
        for j in range(Nx): #space: xs[j]
        
            # Find h(t,x) at x = [xs[j-1], xs[j], xs[j+1]] = [hl, hc, hr]
            hl = hs[i][(j-1) % Nx] 
            hc = hs[i][j % Nx]       
            hr = hs[i][(j+1) % Nx] 
            
            D_center[j] = hr**3 + 3*hc * (hr**2 + hr*hc + hc*hl + hl**2) + 2*hc**3 + hl**3 # * -1, below
            
            if j < Nx-1:
                D_upper[j] = hr**3 + 3*hc*hr * (hc + hr) + hc**3
                
            if j > 0:
                D_lower[j-1] = hl**3 + 3*hc*hl * (hc + hl) + hc**3
        
        # combine as upper, middle, lower diagonals
        D = np.diagflat(-D_center) + np.diagflat(D_lower, -1) + np.diagflat(D_upper, 1)
    
    
        # adjust for periodic boundary...
        
        # -- set top right corner to D_lower with j = 0
        
        hl = hs[i][Nx-1] 
        hc  = hs[i][0]       
        hr = hs[i][1] 
        
        D[0, Nx-1] = hl**3 + 3*hc*hl * (hc + hl) + hc**3
        
        # -- set bottom left corner to D_upper at j = N-1 
        
        hl = hs[i][Nx-2] 
        hc  = hs[i][Nx-1]       
        hr = hs[i][0] 
        
        D[Nx-1, 0] = hr**3 + 3*hc*hr * (hc + hr) + hc**3

        # singular finite difference matrix
        D = D / (8 * dx**2)
        
        # assume sum p'' = 0 
        
        for k in range(Nx):
            D[Nx-1, k] = 1
            
        f_n[i][Nx-1] = 0
        
        # solve for p
        p_n = np.linalg.solve(D, f_n[i])
    
    
        #return D, p_n
    
        
        # plotting ( still inside time loop )
        if (i*figrs) % Nt == 0: 
            
            pp.figure()
            
            pp.plot(xs, p_exact, label="$p$", color='r')
            
            pp.plot(xs, p_n, label="$p_N$", color='b');
            
            #plot h^3 u'
            # hs_cube = [h**3 for h in hs[i]]
            # p_dxs = [p_dx(x) for x in xs]
            # inner = [hs_cube[i] * p_dxs[i] for i in range(Nx)]
            # pp.plot(xs, inner, label="$h^3 p'$", color='black');
            
            #plot f
            #pp.plot(xs, f_n[i], label="f", color='g')
            
            pp.xlabel('x')
            pp.ylabel('p(x)')
            pp.legend()
            
            pp.title("$(h^3p')'= f$ | $p = %s$ | $h=%s$ | $N_x=%d$ | $t = %d$"%(str_p, str_h, Nx, ts[i]))
        
        
        # error at t = ts[i]     
        inf_norms[i] = np.max(np.abs(np.subtract(p_exact,p_n)))
        
        print ("Solved Nx=%d with error %0.5f at t=%d"%(Nx, inf_norms[i], ts[i]))    
    return inf_norms

#------------------------------------------------------------------------------
# Convergence
#------------------------------------------------------------------------------   
def conveg(trials=6, N0=5):
    N = N0
    infNorms = np.zeros(trials)
    dxs = np.zeros(trials)
    dxs_sqr = np.zeros(trials)
    
    for i in range(trials):
        
        infNorms[i] = np.max(solve(N, 1)) # max over time
        dxs[i] = ((xb - xa) / N)
        dxs_sqr[i] = dxs[i]**2
        
        N *= 2
    
    pp.figure(trials+1)
    pp.loglog(dxs, dxs_sqr, color='r', label='$dx^2$')
    pp.loglog(dxs, infNorms, color='b', label='Linf Error')
    
    pp.xlabel('dx')

    pp.legend()
    pp.title('$N_x$=%i to $N_x$=%i with %i trials'%(N0, N/2, trials))


















