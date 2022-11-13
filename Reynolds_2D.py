#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 21 09:43:24 2022

@author: sarahdennis
"""

import numpy as np
import scipy.sparse.linalg as spl
from time import perf_counter
from matplotlib import pyplot as pp
import _graphics as graph

# figs = [0: none, 1: numerical & exact pressure, 2: error]
# BCs = [0: periodic, 1: fixed pressure, 2: fixed x & periodic z]

#------------------------------------------------------------------------------
# Domain
#------------------------------------------------------------------------------
x0, xf = [0, 2*np.pi]
z0, zf = [0, 2*np.pi]
    
def make_grid(Nx, Nz, x0, xf, z0, zf):
    dx = (xf-x0)/(Nx-1)
    dz = (zf-z0)/(Nz-1)
    xs = [x0 + i*dx for i in range(Nx)]
    zs = [z0 + j*dz for j in range(Nz)]
    return [xs, zs]

# Evaluate f(x, z) on grid
# returns f[x][z]
def eval_grid(f, xs, zs):
    f_n = np.zeros((len(xs),len(zs))) #f[t][x][z]
    for i in range(len(xs)):
        for j in range(len(zs)):
            f_n[i][j] = f(xs[i], zs[j])
    return f_n

def unravel(A_flat, Nx, Nz):
    A = np.zeros((Nx, Nz))
    for i in range(Nx):
        for j in range(Nz):
            A[i][j] = A_flat[i*Nz + j]
    return A

# -- height --------------------------------------------------------------------

# option 1: sinusoidal height
## h(X) = h0 + delta * sin(k0 x) * cos(k1 z)

h0, h_delta = [0.2, 0.15] # delta < h0 
h_alpha, h_beta = [1, 1]

h_str = "h(x,z) = %0.2f + %0.2f \sin(%d x) \cos(%d z)"%(h0, h_delta, h_alpha, h_beta) 
def h(x, z):
    return h0 + h_delta * np.sin(h_alpha*x) * np.cos(h_beta*z)

def h_dX(x, z): # = [h_x, h_z]
    return [h_delta * h_alpha * np.cos(h_alpha*x)* np.cos(h_beta*z), -h_delta * h_beta * np.sin(h_alpha*x) * np.sin(h_beta*z)]

# option 2: wedge height (must use BC = 1 or 2)
## h(X) = hf + mx

# h0, hf = [0.1, 0.01] # h0 > hf > 0
# h_m = (2*hf - h0)/(xf - x0)

# h_str = "h(x) = %0.2f + %0.2f x"%(h0, h_m)

# def h(x, z):
#     return hf + h_m*(x-xf)

# def h_dX(x, z): # = [h_t, h_x, h_z]
#     return [h_m, 0]

def plot_height():
    xs, zs = make_grid(graph.Nx, graph.Nz, x0, xf, z0, zf)
    h_grid = eval_grid(h, xs, zs)
    graph.plot_3D(h_grid, xs, zs, "Height: $%s$"%h_str)


#------------------------------------------------------------------------------
# Manf. solution
#------------------------------------------------------------------------------
# (h^3 P')' = f(h)

# Pressure P(xy)

p_alpha, p_beta = [1,1]
def p(x, z):
    return -np.sin(p_alpha*x)*np.cos(p_beta*z)

p_str = "p(x, z) = -\sin(%d x)\cos(%d z)"%(p_alpha, p_beta)

def p_dX(x, z): 
    p_x = -p_alpha * np.cos(p_alpha*x) * np.cos(p_beta*z)
    p_z = p_beta * np.sin(p_alpha*x) * np.sin(p_beta*z)
    return [p_x, p_z]

def p_dXX(x, z):
    p_xx = p_alpha**2 * np.sin(p_alpha*x) * np.cos(p_beta*z)
    p_zz = p_beta**2 * np.sin(p_alpha*x) * np.cos(p_beta*z)
    return [p_xx , p_zz]

def plot_pressure():
    xs, zs = make_grid(graph.Nx, graph.Nz, x0, xf, z0, zf)
    p_grid = eval_grid(p, xs, zs)
    graph.plot_3D(p_grid, xs, zs, "Pressure: $%s$"%p_str)
    
#------------------------------------------------------------------------------
# RHS 
#------------------------------------------------------------------------------
# manufactured solution 
# f = h^3 p'' + (h^3)' p' 
def f(x, z): 
    p_x, p_z = p_dX(x, z)
    p_xx, p_zz = p_dXX(x, z)
    h_x, h_z = h_dX(x, z)
    
    return h(x, z)**3 * (p_xx + p_zz) + 3*h(x, z)**2 * (h_x*p_x + h_z*p_z)

# Reynolds RHS
# f = ht + 0.5 ((h Ux)x + (h Uz)z)
# def f(t, X):
#     h_x, h_z = h_dX(t, X)
#     h_t = h_dt(t,X)
#     h_ = h(t,X)

#     return  0.5*h_*(h_x + h_z)

#------------------------------------------------------------------------------
# numerical solution
#------------------------------------------------------------------------------   

# BCs = [0: periodic, 1: fixed pressure, 2: fixed x & periodic z]
# figs = [0: none, 1: numerical & exact pressure, 2: error]
BCs_txt = ['x & z periodic', 'x & z prescribed', 'x prescribed, z periodic']

def solve(Nx=100, Nz=100, fig=2, BC=1):
    start = perf_counter()
    err = 1
    
    # grid 
    xs, zs = make_grid(Nx, Nz, x0, xf, z0, zf)
    dx = (xf-x0)/(Nx-1)
    dz = (zf-z0)/(Nz-1)
    

    # RHS f(x, z) = f_n[xz]
    f_n = eval_grid(f, xs, zs).ravel()
        
    
    # finite difference matrix D
    D = np.zeros((Nx*Nz, Nx*Nz))
    
    for i in range(Nx): # xs[i]

        # Build rows [i*Nz : (i+1)*Nz]
        
        # Initialize 5 diagonals... 
        P_center = np.zeros(Nz)   #P(i,j)    central diagonal 
        P_west = np.zeros(Nz)     #P(i, j-1) lower diagonal
        P_east = np.zeros(Nz)     #P(i, j+1) upper diagonal 
        P_south = np.zeros(Nz)    #P(i-1, j) left diagonal
        P_north = np.zeros(Nz)    #P(i+1, j) north diagonal
    
        for j in range(Nz): # zs[j]
            
            # Find h(t,X) on grid                  
            h_c = h(xs[i], zs[j])         #h(i,j)    
            h_N = h(xs[(i+1)%Nx], zs[j])  #h(i+1, j)
            h_S = h(xs[(i-1)%Nx],  zs[j]) #h(i-1, j)
            h_E = h(xs[i], zs[(j+1)%Nz])  #h(i, j+1)
            h_W = h(xs[i], zs[(j-1)%Nz])  #h(i, j-1)

            # h^3 = (h^3 + h^3)/2
            #P(i,j) central diagonal
            P_center_x = -(h_N**3 + 2*h_c**3 + h_S**3)/(2* dx**2) #Pij from Px
            P_center_z = -(h_E**3 + 2*h_c**3 + h_W**3)/(2* dz**2) #Pij from Py
            P_center[j] = P_center_x + P_center_z
            
            #P(i+1, j) right diagonal
            P_north[j] = (h_c**3 + h_N**3)/(2*dx**2)
            
            #P(i-1, j) left diagonal
            P_south[j] = (h_c**3 + h_S**3)/(2*dx**2)
            
            #P(i, j-1) lower diagonal
            P_west[j]  = (h_c**3 + h_W**3)/(2*dz**2)
            
            #P(i, j+1) upper diagonal 
            P_east[j]  = (h_c**3 + h_E**3)/(2*dz**2)
            

        #input rows [i*Nx : (i+1)*Nx] into D, adjust for boundary conditions
        
        if BC == 0: #periodic 
        
            # Make three (Nz * Nz) matrices with these diagonals
            D_i_left = np.diagflat(P_south)
            D_i_center = np.diagflat(P_west[1:Nz], -1) + np.diagflat(P_center) + np.diagflat(P_east[0:Nz-1], 1)
            D_i_right = np.diagflat(P_north)
            
            # adjust for periodic boundary in z
            D_i_center[0][Nz-1] = P_west[0]
            D_i_center[Nz-1][0] = P_east[Nz-1]

            # enter (Nz * Nz) blocks into D
            if i == 0:
                D[0:Nz, 0:Nz] = D_i_center
                D[0:Nz, Nz:2*Nz] = D_i_right
                D[0:Nz, (Nx-1)*Nz:Nx*Nz] = D_i_left
            elif i == Nx-1:
                D[(Nx-1)*Nz:Nx*Nz, 0:Nz] = D_i_right
                D[(Nx-1)*Nz:Nx*Nz, (Nx-2)*Nz:(Nx-1)*Nz] = D_i_left
                D[(Nx-1)*Nz:Nx*Nz, (Nx-1)*Nz:Nx*Nz] = D_i_center
            else:
                D[i*Nz:(i+1)*Nz, (i-1)*Nz:i*Nz] = D_i_left
                D[i*Nz:(i+1)*Nz, i*Nz:(i+1)*Nz] = D_i_center
                D[i*Nz:(i+1)*Nz, (i+1)*Nz:(i+2)*Nz] = D_i_right
            
            # closure: assume sum p'' = 0 
            D[Nx*Nz - 1] = 1 # [1 ... 1]
            f_n[Nx*Nz - 1] = 0
        
        elif BC == 1: #fixed in x and z
        
            #set pressure at x0 and xf
            if i == 0 or i == Nx-1: 
            
                P_center = np.ones(Nz)
                
                P_west = np.zeros(Nz)   
                P_east = np.zeros(Nz)     
                P_south = np.zeros(Nz)    
                P_north = np.zeros(Nz)
                
                f_n[i*Nz:(i+1)*Nz] = [p(xs[i], zs[s]) for s in range(0, Nz)]
                
            #set pressure at z0 and zf  
            else: 
                P_center[0] = 1
                P_center[Nz-1] = 1
                
                P_west[0] = 0
                P_west[Nz-1] = 0
                P_east[0] = 0
                P_east[Nz-1] = 0
                P_south[0] = 0
                P_south[Nz-1] = 0
                P_north[0] = 0
                P_north[Nz-1] = 0
                
                f_n[i*Nz] = p(xs[i], zs[0])
                f_n[i*Nz + Nz-1] = p(xs[i], zs[Nz-1])

            # Make three (Nz * Nz) matrices with these diagonals
            D_i_left = np.diagflat(P_south)
            D_i_center = np.diagflat(P_west[1:Nz], -1) + np.diagflat(P_center) + np.diagflat(P_east[0:Nz-1], 1)
            D_i_right = np.diagflat(P_north)
            
            # input rows [i*Nx : (i+1)*Nx] into D
            if i == 0 or i == Nx-1: 
                D[i*Nz:(i+1)*Nz, i*Nz:(i+1)*Nz] = D_i_center
                
            else:
                D[i*Nz:(i+1)*Nz, (i-1)*Nz:i*Nz] = D_i_left
                D[i*Nz:(i+1)*Nz, i*Nz:(i+1)*Nz] = D_i_center
                D[i*Nz:(i+1)*Nz, (i+1)*Nz:(i+2)*Nz] = D_i_right
                
        elif BC == 2: #fixed pressure in x, periodic in z
        
            #fix pressure at x0 and xf
            if i == 0 or i == Nx-1: 
            
                P_center = np.ones(Nz)
                
                P_west = np.zeros(Nz)   
                P_east = np.zeros(Nz)     
                P_south = np.zeros(Nz)    
                P_north = np.zeros(Nz)
                
                f_n[i*Nz:(i+1)*Nz] = [p(xs[i], zs[j]) for j in range(0, Nz)] #p(xa, zs[0:Nz])
           
            # Make three (Nz * Nz) matrices with these diagonals
            D_i_left = np.diagflat(P_south)
            D_i_center = np.diagflat(P_west[1:Nz], -1) + np.diagflat(P_center) + np.diagflat(P_east[0:Nz-1], 1)
            D_i_right = np.diagflat(P_north)
                
            #input rows into D
            if i == 0 or i == Nx-1:
                D[i*Nz:(i+1)*Nz, i*Nz:(i+1)*Nz] = D_i_center
            else: 
                # adjust for periodic boundary in z
                D_i_center[0][Nz-1] = P_west[0]
                D_i_center[Nz-1][0] = P_east[Nz-1]
                
                D[i*Nz:(i+1)*Nz, (i-1)*Nz:i*Nz] = D_i_left
                D[i*Nz:(i+1)*Nz, i*Nz:(i+1)*Nz] = D_i_center
                D[i*Nz:(i+1)*Nz, (i+1)*Nz:(i+2)*Nz] = D_i_right
                
        
    # ------ D is ready now! ----------------------------------
    
    # solve for p_n = [p00, p01, p02, ..., p0N, p10, p11, ..., pNN]
    solve_start = perf_counter()
    p_n = spl.spsolve(D, f_n)
    solve_end = perf_counter()
    
    # unravel p_n for graphing
    p_n_2D = unravel(p_n, Nx, Nz)

    # error  
    p_exact_2D = eval_grid(p, xs, zs)
    err = np.max(np.abs(np.subtract(p_exact_2D,p_n_2D)))
    
    # Plotting
    if fig == 1: # p_n(x, z)
        title = "$%s$ | $%s$ \n $N_x, N_z=%d, %d$ | BC: %s"%(p_str, h_str, Nx, Nz, BCs_txt[BC])
        graph.plot_3D(p_n_2D, xs, zs, title)
    
    elif fig==2: # Error
        title = "Error: $p_n - p_exact$ | $N_x=%d$ , $N_z=%d$ | $dx =%.3f$ $dz = %.3f$ \n | $%s$ | BC: %s"%(Nx, Nz, dx, dz, h_str, BCs_txt[BC])
        graph.plot_3D(p_n_2D-p_exact_2D, xs, zs, title)
       
    print("Solved x:[%.1f,%.1f] z:[%.1f, %.1f] with (Nx=%d, Nz=%d) and BC: %s"%(x0, xf, z0, zf, Nx, Nz, BCs_txt[BC]))    
    print("Max error %0.5f"%(err))
    
    end = perf_counter()
    
    print("Solve time = %0.3f"%(solve_end-solve_start))
    print("Run time = %0.3f"%(end-start))
    #print("solve/run ratio %0.4f"%((solve_end-solve_start)/(end-start)))
    print("\n")
    return err

#------------------------------------------------------------------------------
# Convergence
#------------------------------------------------------------------------------ 
# BCs = [0: periodic, 1: fixed pressure, 2: fixed x & periodic z]
def conveg(trials=5, N0=10, BC=1):
    trial_Nx = N0
    trial_Nz = N0
    
    infNorms = np.zeros(trials)
    
    dxs = np.zeros(trials)
    dxs_sqr = np.zeros(trials)
    
    fig = 0

    for i in range(trials):
        
        if i == trials-1: #make 3D plot at the last iteration
            fig = 2
       
        print('starting trial %d of %d'%(i+1, trials))
            
        infNorms[i] = solve(trial_Nx, trial_Nx, fig, BC)
        
        dx = (xf - x0) / (trial_Nx-1)
        dz = (zf - z0) / (trial_Nz-1)
        dxs[i] = dx + dz
        dxs_sqr[i] = (dx**2) + (dz**2)
        
        trial_Nx *= 2
        trial_Nz *= 2
    
    pp.figure(trials+1)
    pp.loglog(dxs, dxs_sqr, color='r', label='$dx^2 + dz^2$')
    pp.loglog(dxs, dxs, color='g', label='$dx + dz$')
    pp.loglog(dxs, infNorms, color='b', label='$L_\inf$ Error')
    pp.xticks(dxs[::2])
    pp.xlabel('$dx$')

    pp.legend()
    pp.title('$N_x$=%i to $N_x$=%i over %i trials | BC: %s \n $%s$ | $%s$'%(N0, trial_Nx*2**trials, trials, BCs_txt[BC], h_str, p_str))
