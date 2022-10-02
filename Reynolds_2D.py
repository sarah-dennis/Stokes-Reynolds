#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 21 09:43:24 2022

@author: sarahdennis
"""

import numpy as np
import scipy.sparse.linalg as spl
from matplotlib import pyplot as pp
from time import perf_counter

#------------------------------------------------------------------------------
# Domain & Height
#------------------------------------------------------------------------------

# space & time
xa, xb = (0, 4*np.pi)
ya, yb = (0, 2*np.pi)
t0, tf = (0, 2*np.pi)

# height h(t,x) = h(x)
  
h0, delta = (0.2, 0.15) # delta < h0 so height positive
k1, k2 = (2, 5) 

def h(t, X):
    x, y = X
    return h0 + delta * np.sin(k1*x) *np.cos(k2*y)

str_h = "%0.2f + %0.2f \sin(%d x) \cos(%d y)"%(h0, delta, k1, k2) #for graph title


def h_dX(t, X): # = [hx, hy]
    x, y = X
    return [delta * k1 * np.cos(k1*x)* np.cos(k2*y), -delta * k2 * np.sin(k1*x) * np.sin(k2*y)]

def eval_height(ts, xs, ys): 
    Nt = len(ts)
    Nx = len(xs)
    Ny = len(ys)
    hs = np.zeros((Nt, Nx, Ny))
    for k in range(Nt): 
        for i in range(Nx):
            for j in range(Ny):
                hs[k][i][j] = h(ts[k], (xs[i], ys[j]))
    return hs
    

#------------------------------------------------------------------------------
# Manf. solution & RHS
#------------------------------------------------------------------------------
# (h^3 P')' = f(h)

# Pressure P(xy)

alpha, beta = (3,-2)
def p(X):
    x, y = X
    return -np.sin(alpha*x)*np.cos(beta*y)

str_p = "-\sin(%d x)\cos(%d y)"%(alpha, beta) #for graph title

def eval_p(xys):
    Nxy = len(xys)
    p_exact = np.zeros(Nxy)
    for i in range(Nxy):
        p_exact[i] = p(xys[i])
    return p_exact

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
    
    return h(t,X)**3 * (pxx + pyy) + 3*h(t,X)**2 * (hx*px + hy*py)

def eval_f(ts, xys):
    Nt = len(ts)
    Nxy = len(xys)
    f_n = np.zeros((Nt, Nxy)) #f[t][xy]
    for i in range (Nt):
        for j in range(Nxy):
            f_n[i][j] = f(ts[i], xys[j])
    return f_n
    


#------------------------------------------------------------------------------
# numerical solution
#------------------------------------------------------------------------------   
# plot settings
view_theta = 30 # pan view angle up/down
view_phi =60  # pan view angle left-right
# figs= [0: return D p_n, 1: slice (x, y0), 2: slice (x0, y), 3: plot 3D]

def solve(Nx=100, Ny=100, Nt = 1, fig=3):
    start = perf_counter()
    
    # grid 
    dt = (tf - t0)/Nt
    dx = (xb - xa)/Nx
    dy = (yb - ya)/Ny

    ts = [t0 + i*dt for i in range(Nt)]
    xs = [xa + i*dx for i in range(Nx)]
    ys = [ya + i*dy for i in range(Ny)]
    
    # flat space [(x0,y0), (x0, y1), ..., (x0, yN), (x1, y0), ..., (xN, yN)]
    xys = np.zeros((Nx*Ny,2)) 
    for i in range(len(xys)):
        xys[i] = [xs[i//Nx], ys[i%Ny]]  

    # height h(t,X) = h[t][x][y]
    hs = eval_height(ts, xs, ys)
    
    # RHS f(t,X) = f_n[t][xy]
    f_n = eval_f(ts, xys)
        
    # exact solution p(t, X) = p[xys]
    p_exact = eval_p(xys)
    
    # Loop over time t = ts[i]
    inf_norms = np.zeros(Nt) 
    for k in range(Nt): #time: ts[k]
    
        # construct finite difference matrix D
        D = np.zeros((Nx*Ny, Nx*Ny))
        
        for i in range(Nx): # xs[i]
            
            # Build rows [i*Ny : (i+1)*Ny] of D
            
            # Initialize 5 diagonals... 
            P_center = np.zeros(Ny)  #P(i,j)    central diagonal 
            P_west = np.zeros(Ny)    #P(i, j-1) lower diagonal
            P_east = np.ones(Ny)     #P(i, j+1) upper diagonal 
            P_south = np.ones(Ny)    #P(i-1, j) left diagonal
            P_north = np.ones(Ny)    #P(i+1, j) north diagonal
        
            for j in range(Ny): # ys[j]
            
                # Find h(t,X) on grid                  
                h_c = hs[k][i][j]        #h(i,j)    
                h_N = hs[k][(i+1)%Nx][j] #h(i+1, j)
                h_S = hs[k][(i-1)%Nx][j] #h(i-1, j)
                h_E = hs[k][i][(j+1)%Ny] #h(i, j+1)
                h_W = hs[k][i][(j-1)%Ny] #h(i, j-1)
                
                #P(i,j) central diagonal
                P_center_x = -(h_N**3 + 3*h_c*h_N**2 + 3*h_N*h_c**2 + 2*h_c**3 + 3*h_S*h_c**2 + 3*h_c*h_S**2 + h_S**3)/(8* dx**2) #Pij from Px
                P_center_y = -(h_E**3 + 3*h_c*h_E**2 + 3*h_E*h_c**2 + 2*h_c**3 + 3*h_W*h_c**2 + 3*h_c*h_W**2 + h_W**3)/(8* dy**2) #Pij from Py
                P_center[j] = P_center_x + P_center_y
                
                #P(i+1, j) right diagonal
                P_north[j] = (h_c**3 + 3*h_c*h_N**2 + 3*h_N*h_c**2 + h_N**3)/(8*dx**2)
                
                #P(i-1, j) left diagonal
                P_south[j] = (h_c**3 + 3*h_c*h_S**2 + 3*h_S*h_c**2 + h_S**3)/(8*dx**2)
                
                #P(i, j-1) lower diagonal
                P_west[j]  = (h_c**3 + 3*h_c*h_W**2 + 3*h_W*h_c**2 + h_W**3)/(8*dy**2)
                
                #P(i, j+1) upper diagonal 
                P_east[j]  = (h_c**3 + 3*h_c*h_E**2 + 3*h_E*h_c**2 + h_E**3)/(8*dy**2)

            
            # Make three (Ny * Ny) matrices with these diagonals
            D_i_left = np.diagflat(P_south)
            D_i_center = np.diagflat(P_west[1:Ny], -1) + np.diagflat(P_center) + np.diagflat(P_east[0:Ny-1], 1)
            D_i_right = np.diagflat(P_north)
            
            # adjust for periodic boundary in y
            D_i_center[0][Ny-1] = P_west[0]
            D_i_center[Ny-1][0] = P_east[Ny-1]
            
            # input rows [i*Nx : (i+1)*Nx] into D
            # adjust for periodic boundary in x
            
            if i == 0:
                D[0:Ny, 0:Ny] = D_i_center
                D[0:Ny, Ny:2*Ny] = D_i_right
                D[0:Ny, (Nx-1)*Ny:Nx*Ny] = D_i_left
            elif i == Nx-1:
                D[(Nx-1)*Ny:Nx*Ny, 0:Ny] = D_i_right
                D[(Nx-1)*Ny:Nx*Ny, (Nx-2)*Ny:(Nx-1)*Ny] = D_i_left
                D[(Nx-1)*Ny:Nx*Ny, (Nx-1)*Ny:Nx*Ny] = D_i_center
            else:
                D[i*Ny:(i+1)*Ny, (i-1)*Ny:i*Ny] = D_i_left
                D[i*Ny:(i+1)*Ny, i*Ny:(i+1)*Ny] = D_i_center
                D[i*Ny:(i+1)*Ny, (i+1)*Ny:(i+2)*Ny] = D_i_right

        # assume sum p'' = 0 
        D[Nx*Ny - 1] = 1 # [1 ... 1]
        f_n[k][Nx*Ny - 1] = 0
        

        
        # solve for p_n = [p00, p01, p02, ..., p0N, p10, p11, ..., pNN]
        solve_start = perf_counter()
        p_n = spl.spsolve(D, f_n[k]) 
        #p_n = np.linalg.solve(D, f_n[k])
        solve_end = perf_counter()
        
        
        # Plotting at t=ts[0]
        if k == 0: 
        
            if fig == 1: #(x, y0, p)
                y0 = ya
                p_n_1D_x = np.zeros(Nx)
                p_exact_1D_x = np.zeros(Nx)

                for i in range(Nx):
                    p_n_1D_x[i] = p_n[i*Nx + y0]
                    p_exact_1D_x[i] = p_exact[i*Nx+y0]
                    
                pp.figure()
                pp.plot(xs, p_n_1D_x, label="$p_n(y_0)$")
                pp.plot(xs, p_exact_1D_x, label="$p_n(y_0)$")
                pp.title("$p= %s$ | $h=%s$ | $N_x=%d$"%(str_p, str_h, Nx))
                pp.xlabel('x')
                pp.ylabel('p')
                pp.legend()

            elif fig == 2: #(x0, y, p)
                x0 = xa
                p_n_1D_y = p_n[x0*Nx:x0*Nx + Ny]
                p_exact_1D_y = p_exact[x0*Nx:x0*Nx + Ny]
                    
                pp.figure()
                pp.plot(ys, p_n_1D_y, label="$p_n(x_0)$")
                pp.plot(ys, p_exact_1D_y, label="$p(x_0)$")
                pp.title("$p= %s$ | $h=%s$ | $N_y=%d$"%(str_p, str_h, Ny))
                pp.xlabel('y')
                pp.ylabel('p')
                pp.legend()
        
            elif fig == 3: #(x, y, p)
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
    
        # error at t = ts[k]     
        inf_norms[k] = np.max(np.abs(np.subtract(p_exact,p_n)))
        
        print("Solved Nx=%d, Ny=%d"%(Nx, Ny))    
        print("Error %0.5f at t=%0.2f"%(inf_norms[k], ts[0]))
        
        end = perf_counter()
        
        print("Total run time = %0.3f"%(end-start))
        print("Solve time = %0.3f"%(solve_end-solve_start))
        #print("solve/run ratio %0.4f"%((solve_end-solve_start)/(end-start)))
   
    return inf_norms

#------------------------------------------------------------------------------
# Convergence
#------------------------------------------------------------------------------   
def conveg(trials=6, N0=5):
    Nx = N0
    Ny = N0
    Nt = 1
    
    infNorms = np.zeros(trials)
    dxs = np.zeros(trials)
    dxs_sqr = np.zeros(trials)
    
    for i in range(trials):
        
        infNorms[i] = np.max(solve(Nx, Ny, Nt, 0)) # max over time
        dxs[i] = ((xb - xa) / Nx)
        dxs_sqr[i] = dxs[i]**2
        
        Nx += 50
        Ny += 50
    
    pp.figure(trials+1)
    pp.loglog(dxs, dxs_sqr, color='r', label='$dx^2$')
    pp.loglog(dxs, infNorms, color='b', label='Linf Error')
    
    pp.xlabel('dx')

    pp.legend()
    pp.title('$N_x$=%i to $N_x$=%i with %i trials'%(N0, Nx-50, trials))

