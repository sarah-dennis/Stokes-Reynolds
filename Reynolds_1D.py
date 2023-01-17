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
# Fluid constants
#------------------------------------------------------------------------------
# viscosity
eta = 1 

#------------------------------------------------------------------------------
# Domain & Height
#------------------------------------------------------------------------------
# width
xa, xb = (0, 2)

# lower surface velocity
U = 2

# Height

# -- option 1: Sinusoidal height 
# h0, delta, k = (1, 0.3, 6) # delta < h0 so height positive
# h_str = "h(x) = %0.1f + %0.1f \cos(%d x)"%(h0, delta, k) #for graph title

# def h(x):
#     return h0 + delta * np.sin(k*x)

# def h_dx(x):
#     return delta * k * np.cos(k*x)

#p0 = 0
#pf = 0

# -- option 2: Wedge height 
# h0, hf = (1, 0.8) #h0 = inital height, hf: final height
# delta = (hf - h0)/(xb - xa) # appropriate slope
# h_str = "h(x) = %0.1f + %0.2f (x - %d)"%(h0, delta, xa) #for graph title
   
# def h(x):
#     return delta * (x - xa) + h0

# def h_dx(x):
#     return delta

          
# p0 = 0
# pf = 0
   
# # -- option 3: Rayleigh step height
# h1, h2 = (1, 0.1)
# xm = (xb-xa)/2
# l1 = xm - xa
# l2 = xb - xm

# h_str = "h(x_0:x_m) = %0.1f, h(x_M:x_N) = %0.1f"%(h1, h2) 

# def h(x):
#     if x <= xm:
#         return h1
#     else:
#         return h2
    
# def h_dx(x):
#     return 0

# p0 = 0

# # max presure for Rayleigh step
# def p_max(h1, h2, l1, l2, p0):
#     return (6*eta*U*(h1 - h2)*l1) / (h1**3 + h2**3*(l1/l2)) + p0

# def p_min(h1, h2, l1, l2, p0):
#     return (-6*eta*U*(h1 - h2)*l1) / (h1**3 + h2**3*(l1/l2)) + p0

# pm = p_max(h1, h2, l1, l2, p0)
# pf = p_min(h1, h2, l1, l2, pm)

# -- option 4: 3-step height
h1, h2, h3 = (0.5, 0.3, 0.4)
xm1 = (xb-xa)/3
xm2 = 2*(xb-xa)/3
l1 = xm1 - xa
l2 = xm2 - xm1
l3 = xb - xm2

h_str = "h(x_0:x_1) = %0.1f, h(x_1:x_2) = %0.1f, h(x_3:x_3) = %0.1f"%(h1, h2, h3) 

def h(x):
    if x <= xm1:
        return h1
    elif x <= xm2:
        return h2
    else: 
        return h3
    
def h_dx(x): 
    return 0 # for manf_rhs()

p0 = 0 #p0
pf = -1 #p3

# max presure for Rayleigh step

p2 = (h1**3*h2**3*l3**2 + h1**3 *h3**3*l2*l3 - h1**6*l2*l3 -h1**3*h2**3*l1*l3)*p0 - (h1**3*h3**3*l1*l2 - h2**3*h3**3*l1**2 - h1**6*l2*l3 + h1**3*h2**3*l1*l3)*pf -6*U*eta*(h1**3*h3*l1*l2*l3 - h1**4*l1*l2*l3 - h1**3*h2*l1*l3*(h2**3*l3+h3**3*l2) + h1**4*l1*l3*(h2**3*l3+h3**3*l2) + h1**3*h2*l3**2*(h1**3*l2+h2**3*l1) -h1**3*h3*l3**2*(h1**3*l2+h2**3*l1)+h2**3*l1**2*l3*(h3-h1))/(h1**3*h2**3*l3**2+h1**3*h3**3*l2*l3-l1*l2*h1**3*h3**3-l1**2*h2**2*h3**3)


p1 = (6*U*eta*(l3*(h2-h3)*(h1**3*l2+h2**3*l1)-l1*(h2-h1)*(h2**3*l3+h3**3*l2))-(h2**3*l3+h3**3*l2)*p0 + (h2**3*l3 + h3**3*l2)*p2 + (h1**3*l2 + h2**3*l1)*pf)/(h1**3*l2+h2**3*l1)

print("extrema: %.2f, %.2f, %.2f, %.2f"%(p0, p1, p2, pf))
#------------------------------------------------------------------------------
# RHS: Manf & Exact solution
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
    return f_n

#------------------------------------------------------------------------------
# RHS: Reynolds equation 
#------------------------------------------------------------------------------
# # RHS = 6 eta U hx

def discr_hx(Nx, dx, h, xs, BC):
    D_lower = -1*np.ones(Nx)
    D_upper = np.ones(Nx)
    D = np.diagflat(D_lower[1:Nx], -1) + np.diagflat(D_upper[0:Nx-1], 1)
    
    if BC == 0: #periodic 
       
        D[0][Nx-1] = -1
        D[Nx-1][0] = 1
        
    if BC == 1: #prescribed 
        None
        # boundary heights are not used
        
    D = D/(2*dx)
    
    hs = np.zeros(Nx) #h[x]
    for i in range(Nx):
        hs[i] = h(xs[i])
        
    hs_dx = D@hs 

    #graph.plot_2D_multi([hs_dx, hs], xs, "height", ["hx", "h"])
    return hs_dx

p_str = "p(x) \; unknown"

#------------------------------------------------------------------------------
# numerical solution
#------------------------------------------------------------------------------   
#figrs = [Exact and numerical Pressure, Error]

#BCs = [0: periodic, 1: fixed]

def solve(Nx=100, figrs=1, BC=1, exact_sol=False):
    
    if BC==0: #periodic
        dx = (xb - xa)/(Nx)
        
    elif BC==1: #fixed
        dx = (xb - xa)/(Nx-1)

    xs = [xa + i*dx for i in range(Nx)]

    # construct RHS on grid
    if exact_sol: 
        inf_norm_err = 0
        f_n = manf_rhs(Nx, xs)
    
        if BC == 1: # set boundary pressures to exact solution
            pa = p(xs[0])
            pb = p(xs[-1])
        
        
    else: #Reynolds RHS
        hx_n = discr_hx(Nx, dx, h, xs, BC)
        f_n = 6 * eta * U * hx_n
        
        if BC == 1: # set boundary pressures to something!
            pa = p0
            pb = pf


    # initilise diagonals of differnce matrix
    D_lower = np.ones(Nx)
    D_center = np.ones(Nx)
    D_upper = np.ones(Nx)
    
    for i in range(Nx): #space: xs[i]
    
        # Find h(t,x) at x = [xs[i-1], xs[i], xs[i+1]] = [hl, hc, hr]
        hl = h(xs[(i-1) % Nx])
        hc = h(xs[i % Nx])   
        hr = h(xs[(i+1) % Nx] )

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
    
    # Plotting and error 
    
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
            label="Error"
            graph.plot_2D(p_exact-p_n, xs, title, label)
        
        return inf_norm_err
    
    else: 
        if figrs == 1:
            title = "Numerical Pressure | $N_x=%d$, $dx = %.2f$ \n $%s$ "%(Nx, dx, h_str)
            labels = ["$p(x)$", "$h(x)$"]
            hs = [h(x) for x in xs]
            graph.plot_2D_multi([p_n, hs] , xs, title, labels)
        
        #max pressure check for step height
        print("Max numerical pressure: %.3f"%np.max(p_n))
        print("Min numerical pressure: %.3f"%np.min(p_n))
        
        #return p_n

#------------------------------------------------------------------------------
# Convergence
#------------------------------------------------------------------------------   
def conveg(trials=10, N0=5, BC=1, figrs=0):
    N = N0
    infNorms = np.zeros(trials)
    dxs = np.zeros(trials)
    dxs_sqr = np.zeros(trials)
    
    for i in range(trials):
        
        infNorms[i] = solve(N, figrs, BC, exact_sol=True)
        
        if BC==0: #periodic
            
            dx = (xb - xa)/(N)
            
        elif BC==1: #fixed
        
            dx = (xb - xa)/(N-1)
        
        dxs[i] = dx
        dxs_sqr[i] = dxs[i]**2
        
        N *= 2
    
    pp.figure(trials+1)
    pp.loglog(dxs, dxs_sqr, color='r', label='$dx^2$')
    pp.loglog(dxs, infNorms, color='b', label='Linf Error')
    
    pp.xlabel('dx')

    pp.legend()
    pp.title('$N_x$=%i to $N_x$=%i with %i trials'%(N0, N/2, trials))
    
    