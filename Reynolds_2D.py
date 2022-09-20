#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 21 09:43:24 2022

some changes

@author: sarahdennis
"""

import numpy as np
from scipy.sparse import linalg as sp
from matplotlib import pyplot as pp
from time import perf_counter

#------------------------------------------------------------------------------
# Domain & Height
#------------------------------------------------------------------------------
# width
xa, xb = (0, 2*np.pi)
ya, yb = (0, 2*np.pi)

#time
t0, tf = (0, 2*np.pi)

# height h(t,x) = h(x)  
h0, delta, k = (0.5, 0.3, 2) # delta < h0 so height positive

str_h = "%0.1f + %0.1f \sin(%d xy)"%(h0, delta, k) #for graph title

def h(t, X):
    x, y = X
    return h0 + delta * np.sin(k*x*y) + x*y

def h_dX(t, X):
    x, y = X
    return [delta*k * y * np.cos(k*x*y) + y, delta*k * x * np.cos(k*x*y) + x]
                
#------------------------------------------------------------------------------
# Manf. solution & RHS
#------------------------------------------------------------------------------
# (h^3 p')' = f(h)

# set "exact" p(X)

alpha = 3
def p(X):
    x, y = X
    return -np.sin(alpha*x*y)

str_p = "-\sin(%d xy)"%(alpha) #for graph title

def p_dX(X): # = [p_x, p_y]
    x, y = X
    return [-(alpha * y) * np.cos(alpha*x*y), -(alpha * x) * np.cos(alpha*x*y)]

def p_dXX(X): # = [p_xx, p_yy]
    x, y = X
    return [(alpha * y)**2 * np.sin(alpha*x*y), (alpha * x)**2 * np.sin(alpha*x*y)]

# calculate RHS f(x,t)

def f(t,X): # = h^3 p'' + (h^3)' p' 
    px, py = p_dX(X)
    pxx, pyy = p_dXX(X)
    hx, hy = h_dX(t, X)
    
    return (h(t,X)**3) * (pxx + pyy) + 3*h(t,X)**2 * (hx*px + hy*py)


#------------------------------------------------------------------------------
# numerical solution
#------------------------------------------------------------------------------   

def solve(NX=100, Nt = 1, figs=1):
    start = perf_counter()
    
    # grid sizing
    Ny = NX
    Nx = NX
    dx = (xb - xa)/Nx
    dy = (yb - ya)/Ny
    dt = (tf - t0)/Nt
    
    # time t
    ts = [t0 + i*dt for i in range(Nt)]
    
    # space X = (x, y) 
    xs = [xa + i*dx for i in range(Nx)]
    ys = [ya + i*dy for i in range(Ny)]
    
    # ravel X -> xys = [(x0,y0), (x0, y1), ..., (x0, yN), (x1, y0), ..., (xN, yN)]
    xys = np.zeros((Nx*Ny,2)) 
    for i in range(len(xys)):
        
        xys[i] = [xs[i//Nx], ys[i%Ny]]

    #height h(t,X)
    hs = np.zeros((Nt, Nx, Ny)) #h[t][x][y]
    for k in range(Nt): 
        for i in range(Nx):
            for j in range(Ny):
                hs[k][i][j] = h(ts[k], (xs[i], ys[j]))
    
    #RHS f(t,X)
    f_n = np.zeros((Nt, Nx*Ny)) #f[t][xy]
    for i in range (Nt):
        for j in range(Nx*Ny):
            f_n[i][j] = f(ts[i], xys[j])
        
    #exact p(t, X)
    p_exact = np.zeros(Nx*Ny)
    for i in range(Nx*Ny):
        p_exact[i] = p(xys[i])
         
        
    # Loop over time t = ts[i]
    
    inf_norms = np.zeros(Nt) 
    
    # construct finite difference matrix
    for k in range(Nt): #time: ts[k]
    
        D = np.zeros((Nx*Ny, Nx*Ny))
        
        for i in range(Nx): # xs[i]
            
            # Prepare to build rows [i*Ny : (i+1)*Ny] of D
            
            # Initialize 5 diagonals... 
            
            # for columns [i*Ny, (i+1)*Ny]
            P_center = np.zeros(Ny) #P(i,j)    central diagonal 
            P_north = np.zeros(Ny)  #P(i, j-1) lower diagonal
            P_south = np.ones(Ny)   #P(i, j+1) upper diagonal 
            
            # for columns [(i-1)*Ny, i*Ny]
            P_east = np.ones(Ny)    #P(i+1, j) right diagonal
            
            # for columns [i*Ny, (i+1)*Ny]
            P_west = np.ones(Ny)    #P(i-1, j) left diagonal
        
            for j in range(Ny): # ys[j]
            
                # Find h(t,X) on grid                  
                h_c = hs[k][i%Nx][j%Ny]     #h(i,j)    
                h_N = hs[k][i%Nx][(j-1)%Ny] #h(i, j-1)
                h_S = hs[k][i%Nx][(j+1)%Ny] #h(i, j+1)
                h_E = hs[k][(i+1)%Nx][j%Ny] #h(i+1, j)
                h_W = hs[k][(i-1)%Nx][j%Ny] #h(i-1, j)
                
                # build diagonals
                
                #P(i,j) central diagonal
                P_center_x = (h_E**3 + 3*h_E**2*h_c + 3*h_E*h_c**2 + 2*h_c**3 + 3*h_c**2*h_W + 3*h_c*h_W**2 + h_W**3) #Pij from Px
                P_center_y = (h_S**3 + 3*h_S**2*h_c + 3*h_S*h_c**2 + 2*h_c**3 + 3*h_c**2*h_N + 3*h_c*h_N**2 + h_N**3) #Pij from Py
                P_center[j] = P_center_x/(8*dx**2) + P_center_y/(8*dy**2) # * -1 below 
                
                #P(i, j-1) lower diagonal
                P_north[(j-1)%Ny] = (h_c**3 + 3*h_N**2*h_c + 3*h_N*h_c**2 + h_N**3)/(8*dy**2)
                
                #P(i, j+1) upper diagonal
                P_south[(j+1)%Ny] = (h_c**3 + 3*h_S**2*h_c + 3*h_S*h_c**2 + h_S**3)/(8*dy**2)
                
                #P(i-1, j) left diagonal
                P_west[j] = (h_c**3 + 3*h_W**2*h_c + 3*h_W*h_c**2 + h_W**2)/(8*dx**2)
                
                #P(i+1, j) right diagonal
                P_east[j] = (h_c**3 + 3*h_E**2*h_c + 3*h_E*h_c**2 + h_E**2)/(8*dx**2)

            
            # Make three (Ny * Ny) matrices with these diagonals
            D_i_left = np.diagflat(P_west)
            D_i_center = np.diagflat(P_north[0:Ny-1], -1) + np.diagflat(-P_center) + np.diagflat(P_south[0:Ny-1], 1)
            D_i_right = np.diagflat(P_east)
            
            # adjust for periodic boundary in y
            D_i_center[0][Ny-1] = P_north[Ny-1]
            D_i_center[Ny-1][0] = P_south[Ny-1]
            
            # input rows [i*Nx : (i+1)*Nx] into D
            # adjust for periodic boundary in x
            if i == 0:
                D[0:Ny, 0:Ny] = D_i_center
                D[0:Ny, Ny:2*Ny] = D_i_right
                D[0:Ny, Nx*Ny-Ny:Nx*Ny] = D_i_left
            elif i == Nx-1:
                D[(Nx-1)*Ny:Nx*Ny, 0:Ny] = D_i_right 
                D[(Nx-1)*Ny:Nx*Ny, (Nx-2)*Ny:(Nx-1)*Ny] = D_i_left
                D[(Nx-1)*Ny:Nx*Ny, (Nx-1)*Ny:Nx*Ny] = D_i_center
            else:
                D[i*Ny:(i+1)*Ny, (i-1)*Ny:i*Ny] = D_i_left
                D[i*Ny:(i+1)*Ny, i*Ny:(i+1)*Ny] = D_i_center
                D[i*Ny:(i+1)*Ny, (i+1)*Ny:(i+2)*Ny] = D_i_right

        # assume sum p'' = 0
        D[Nx*Ny-1] = 1 # = [1 ... 1]
        f_n[k][Nx*Ny - 1] = 0
        
        solve_start = perf_counter()
        # solve for p
        p_n = sp.spsolve(D, f_n[k]) # = [p00, p01, p02, ..., p0N, p10, p11, ..., pNN]
        solve_end = perf_counter()
        
        # plotting ( still inside time loop )
        if k == 0 and figs: #plot at time t = ts[k]
            for i in range(Nx):
                if i == 0 or i == Nx//3 or i == 2*Nx//3 or i == Nx-1: #plot slice at x = xs[i]
                    pp.figure()
                    
                    pp.plot(ys, p_exact[i * Nx: (i+1) * Nx], label="$p$", color='r')
                    
                    pp.plot(ys, p_n[i * Nx: (i+1) * Nx], label="$p_N$", color='b');
                    
                    pp.xlabel('y')
                    pp.ylabel('$p(t = %d, x = %.3f)$'%(ts[k], xs[i]))
                    pp.legend()
                    
                    pp.title("$(h^3p')'= f$ | $p = %s$ | $h=%s$ | $N_x=%d$"%(str_p, str_h, Nx))
            
        
        # error at t = ts[i]     
        inf_norms[k] = np.max(np.abs(np.subtract(p_exact,p_n)))
        
        print ("Solved Nx=%d with error %0.5f at t=%0.2f"%(Nx, inf_norms[k], ts[k]))    
        end= perf_counter()
        
        print("Solve time = %0.3f"%(solve_end-solve_start))
        print("run time = %0.3f"%(end-start))
        print("solve/run ratio %0.4f"%((solve_end-solve_start)/(end-start)))
   
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
        
        infNorms[i] = np.max(solve(N, 1, 0)) # max over time
        dxs[i] = ((xb - xa) / N)
        dxs_sqr[i] = dxs[i]**2
        
        N += 50
    
    pp.figure(trials+1)
    pp.loglog(dxs, dxs_sqr, color='r', label='$dx^2$')
    pp.loglog(dxs, infNorms, color='b', label='Linf Error')
    
    pp.xlabel('dx')

    pp.legend()
    pp.title('$N_x$=%i to $N_x$=%i with %i trials'%(N0, N/2, trials))


















