#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 21 09:43:24 2022

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
xa, xb = (0, np.pi)
ya, yb = (0, 0.5*np.pi)

#time   
t0, tf = (0, 2*np.pi)

# height h(t,x) = h(x)  
h0, delta, k1, k2 = (0.5, 0.1, 4, 1) # delta < h0 so height positive

str_h = "%0.1f + %0.1f \sin(%d x) \cos(%d y)"%(h0, delta, k1, k2) #for graph title

def h(t, X):
    x, y = X
    return h0 + delta * np.sin(k1*x) * np.cos(k2*y)

def h_dX(t, X): # = [hx, hy]
    x, y = X
    return [delta * k1 * np.cos(k1*x)* np.cos(k2*y), -delta * k2 * np.sin(k1*x) * np.sin(k2*y)]
                
#------------------------------------------------------------------------------
# Manf. solution & RHS
#------------------------------------------------------------------------------
# (h^3 p')' = f(h)

# set "exact" p(X)

alpha, beta = (2, 1)
def p(X):
    x, y = X
    return -np.sin(alpha*x)*np.cos(beta*y)

str_p = "-\sin(%d x)\cos(%d y)"%(alpha, beta) #for graph title

def p_dX(X): # = [p_x, p_y]
    x, y = X
    return [-alpha * np.cos(alpha*x) * np.cos(beta*y), beta * np.sin(alpha*x) * np.sin(beta*y)]

def p_dXX(X): # = [p_xx, p_yy]
    x, y = X
    return [alpha**2 * np.sin(alpha*x) * np.cos(beta*y) , beta**2 * np.sin(alpha*x) * np.cos(beta*y)]

# calculate RHS f(x,t)

def f(t,X): # = h^3 p'' + (h^3)' p' 
    px, py = p_dX(X)
    pxx, pyy = p_dXX(X)
    hx, hy = h_dX(t, X)
    
    return (h(t,X)**3) * (pxx + pyy) + 3*h(t,X)**2 * (hx*px + hy*py)


#------------------------------------------------------------------------------
# numerical solution
#------------------------------------------------------------------------------   
# plot settings
view_theta = 10 # pan view angle up/down
view_phi = 20  # pan view angle left-right
figs = True

def solve(NX=100, Nt = 1):
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

    xys = np.zeros((Nx*Ny,2)) #[(x0,y0), (x0, y1), ..., (x0, yN), (x1, y0), ..., (xN, yN)]
    for i in range(len(xys)):
        xys[i] = [xs[i//Nx], ys[i%Ny]]

    # height h(t,X)
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
                h_c = hs[k][i][j]     #h(i,j)    
                h_N = hs[k][i][(j-1)%Ny] #h(i, j-1)
                h_S = hs[k][i][(j+1)%Ny] #h(i, j+1)
                h_E = hs[k][(i+1)%Nx][j] #h(i+1, j)
                h_W = hs[k][(i-1)%Nx][j] #h(i-1, j)
                
                # build diagonals
                
                #P(i,j) central diagonal
                P_center_x = (h_E**3 + 3*h_c*h_E**2 + 3*h_E*h_c**2 + 2*h_c**3 + 3*h_W*h_c**2 + 3*h_c*h_W**2 + h_W**3) #Pij from Px
                P_center_y = (h_S**3 + 3*h_c*h_S**2 + 3*h_S*h_c**2 + 2*h_c**3 + 3*h_N*h_c**2 + 3*h_c*h_N**2 + h_N**3) #Pij from Py
                P_center[j] = -P_center_x/(8*dx**2) - P_center_y/(8*dy**2)
                
                #P(i, j-1) lower diagonal
                P_north[(j-1)%Ny] = (h_c**3 + 3*h_c*h_N**2 + 3*h_N*h_c**2 + h_N**3)/(8*dy**2)
                
                #P(i, j+1) upper diagonal
                P_south[(j+1)%Ny] = (h_c**3 + 3*h_c*h_S**2 + 3*h_S*h_c**2 + h_S**3)/(8*dy**2)
                
                #P(i-1, j) left diagonal
                P_west[j] = (h_c**3 + 3*h_W**2*h_c + 3*h_W*h_c**2 + h_W**2)/(8*dx**2)
                
                #P(i+1, j) right diagonal
                P_east[j] = (h_c**3 + 3*h_E**2*h_c + 3*h_E*h_c**2 + h_E**2)/(8*dx**2)

            
            # Make three (Ny * Ny) matrices with these diagonals
            D_i_left = np.diagflat(P_west)
            D_i_center = np.diagflat(P_north[0:Ny-1], -1) + np.diagflat(P_center) + np.diagflat(P_south[1:Ny], 1)
            D_i_right = np.diagflat(P_east)
            
            # adjust for periodic boundary in y
            D_i_center[0][Ny-1] = P_north[Ny-1] #
            D_i_center[Ny-1][0] = P_south[0]
            
            # input rows [i*Nx : (i+1)*Nx] into D
            # adjust for periodic boundary in x
            if i == 0:
                D[0:Ny, 0:Ny] = D_i_center
                D[0:Ny, Ny:2*Ny] = D_i_right
                D[0:Ny, (Nx-1)*Ny:Nx*Ny] = D_i_left
            elif i == Nx-1:
                D[(Nx-1)*Ny:Nx*Ny, 0:Ny] = D_i_right
                D[(Nx-1)*Ny:Nx*Ny, (i-1)*Ny:(Nx-1)*Ny] = D_i_left
                D[(Nx-1)*Ny:Nx*Ny, i*Ny:Nx*Ny] = D_i_center
            else:
                D[i*Ny:(i+1)*Ny, (i-1)*Ny:i*Ny] = D_i_left
                D[i*Ny:(i+1)*Ny, i*Ny:(i+1)*Ny] = D_i_center
                D[i*Ny:(i+1)*Ny, (i+1)*Ny:(i+2)*Ny] = D_i_right

        # assume sum p'' = 0
        D[Nx*Ny-1] = 1 # [1 ... 1]
        f_n[k][Nx*Ny - 1] = 0
        
        #return D
        
        # solve for p
        solve_start = perf_counter()
        p_n = sp.spsolve(D, f_n[k]) # [p00, p01, p02, ..., p0N, p10, p11, ..., pNN]
        solve_end = perf_counter()
        
        # plotting
        if k == 0 and figs: #plot at time t = ts[k]
                        
            p_n_2D = np.zeros((Nx, Ny))
            p_exact_2D = np.zeros((Nx, Ny))
            
            for i in range(Nx):
                for j in range(Ny):
                    p_n_2D[i,j] = p_n[i*Nx + j]
                    p_exact_2D[i,j] = p_exact[i*Nx + j]
                    
            X, Y = np.meshgrid(xs, ys)        
                    
            pp.figure()
            ax = pp.axes(projection='3d')
            ax.plot_surface(X, Y, p_n_2D, label="$p_N$", rstride=1, cstride=1,cmap='viridis')
            pp.title("$p= %s$ | $h=%s$ | $N_x=%d$"%(str_p, str_h, Nx))
            pp.xlabel('x')
            pp.ylabel('y')
            ax.view_init(view_theta, view_phi)
            
            pp.figure()
            ax = pp.axes(projection='3d')
            ax.plot_surface(X, Y, p_exact_2D, label="$p_exact$", rstride=1, cstride=1, cmap='viridis')
            pp.title("$p = %s$ | $h=%s$ | exact"%(str_p, str_h))
            pp.xlabel('x')
            pp.ylabel('y')
            ax.view_init(view_theta, view_phi)
    
        # error at t = ts[i]     
        inf_norms[k] = np.max(np.abs(np.subtract(p_exact,p_n)))
        
        print ("Solved Nx=%d with error %0.5f at t=%0.2f"%(Nx, inf_norms[k], ts[k]))    
        
        end = perf_counter()
        
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


















