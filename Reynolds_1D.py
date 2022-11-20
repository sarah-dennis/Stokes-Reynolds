#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 21 09:43:24 2022

@author: sarahdennis
"""

import numpy as np
from matplotlib import pyplot as pp
import _graphics as graph
#------------------------------------------------------------------------------
# Domain & Height
#------------------------------------------------------------------------------

# width
xa, xb = (0, 2*np.pi)

# Height

# -- option 1: Sinusoidal height 
# h0, delta, k = (1, 0.3, 6) # delta < h0 so height positive
# h_str = "h(x) = %0.1f + %0.1f \cos(%d x)"%(h0, delta, k) #for graph title

# def h(x):
#     return h0 + delta * np.sin(k*x)

# def h_dx(x):
#     return delta * k * np.cos(k*x)

# -- option 2: Wedge height 
h0, hf = (1, 0.01) #h0 = inital height, hf: final height
delta = (hf - h0)/(xb - xa) # appropriate slope

def h(x):
    return delta * (x - xa) + h0

def h_dx(x):
    return delta

h_str = "h(x) = %0.1f + %0.2f (x - %d)"%(h0, delta, xa) #for graph title
                
# -- option 3: Rayleigh step height
# h0, h1 = (1, 0.5)
# xm = (xb-xa)/2
# l0 = xm - xa
# l1 = xb - xm

# h_str = "h(x_0:x_m) = %0.1f, h(x_M:x_N) = %0.1f"%(h0, h1) 

# def h(x):
#     if x <= xm:
#         return h0
#     else:
#         return h1
    
# def h_dx(x):
#     return 0

#------------------------------------------------------------------------------
# Manf. solution
#------------------------------------------------------------------------------
# # RHS = (h^3 p')'

alpha = 1
def p(x):
    return -np.sin(alpha*x)

p_str = "p(x) = -\sin(%dx)"%(alpha) #for graph title

def p_dx(x):
    return -alpha * np.cos(alpha*x)

def p_dxx(x): 
    return alpha**2 * np.sin(alpha*x)

def manf_rhs(Nx, xs):
    f_n = np.zeros(Nx) #f[x]
    for i in range(Nx):
        f_n[i] = (h(xs[i]) ** 3) * p_dxx(xs[i]) + 3 * h(xs[i])**2 * h_dx(xs[i]) * p_dx(xs[i])

# # boundary pressures set with exact solution
# pa = p(xa)
# pb = p(xb)


#------------------------------------------------------------------------------
# RHS: Reynolds equation 
#------------------------------------------------------------------------------
# (h^3 p')' = 6 eta U hx

eta = 1 # viscosity
U = 1 # lower surface velocity

def discr_hx(Nx, dx, h, xs, BC):
    D_center = -2*np.ones(Nx)
    D_lower = np.ones(Nx)
    D_upper = np.ones(Nx)
    D = np.diagflat(D_center) + np.diagflat(D_lower[1:Nx], -1) + np.diagflat(D_upper[0:Nx-1], 1)
    
    if BC == 0: #periodic 
        D[0][Nx-1] = 1
        D[Nx-1][0] = 1
        
    #if BC == 1: #prescribed 
        # boundary heights are not used
        
    D = D/dx**2
    
    hs = np.zeros(Nx) #h[x]
    for i in range(Nx):
        hs[i] = h(xs[i])
        
    hx_n = np.linalg.solve(D, hs) 

    return hx_n

p_str = "p(x) \; unknown"

# boundary pressures
pa=1
pb=1

#------------------------------------------------------------------------------
# numerical solution
#------------------------------------------------------------------------------   
#figrs = [Exact and numerical Pressure, Error]
#BCs = [periodic, fixed]

def solve(Nx=100, figrs=1, BC=1, exact_sol=False):
    
    if BC==0:
        
        dx = (xb - xa)/(Nx)

        xs = [xa + i*dx for i in range(Nx)]
        
    elif BC==1:
        dx = (xb - xa)/(Nx-1)

        xs = [xa + i*dx for i in range(Nx)]

    # construct RHS on grid
    if exact_sol:
        inf_norm_err = 0
        f_n = manf_rhs(xs)
        
    else:
        hx_n = discr_hx(Nx, dx, h, xs, BC)
        f_n = 6 * eta * U * hx_n


    # initilise diagonals of differnce matrix
    D_lower = np.ones(Nx)
    D_center = np.ones(Nx)
    D_upper = np.ones(Nx)
    
    for i in range(Nx): #space: xs[i]
    
        # Find h(t,x) at x = [xs[i-1], xs[i], xs[i+1]] = [hl, hc, hr]
        hl = h(xs[(i-1) % Nx])
        hc = h(xs[i % Nx])   
        hr = h(xs[(i+1) % Nx] )
        
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
        f_n[0] = pa
        
        # -- set bottom row D to [ ... , 0, 1] and f[Nx-1] = p_outlet
        D[Nx-1,Nx-1] = 1
        D[Nx-1,Nx-2] = 0
        f_n[Nx-1] = pb

    # solve for p
    p_n = np.linalg.solve(D, f_n)
    #return p_n
        

    if exact_sol:
        p_exact = np.zeros(Nx)
        for i in range(Nx):
            p_exact[i] = p(xs[i])
        
        inf_norm_err = np.max(np.abs(np.subtract(p_exact,p_n)))
        print ("Solved Nx=%d with error %0.5f"%(Nx, inf_norm_err))  
        
        if figrs == 1:
            title = "Exact vs. Numerical Pressure \n $%s$ | $%s$"%(p_str,h_str)
            labels = ["exact", "numerical"]
            
            graph.plot_2D_multi([p_exact, p_n], xs, title, labels)
            
        elif figrs == 2:
            title = "Model Error | $N_x=%d$, $dx = %.2f$ \n $%s$ | $%s$"%(Nx, dx, p_str,h_str)
            graph.plot_2D(p_exact-p_n, xs, title)
        
        return inf_norm_err
    
    else: 
        if figrs == 1:
            title = "Numerical Pressure for $%s$"%(h_str)
            graph.plot_2D(p_n, xs, title)
        
            

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



