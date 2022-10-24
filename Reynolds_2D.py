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
# time
t0, tf = [0, 0]

# space
xa, xb = [0, 4*np.pi]
za, zb = [0, 2*np.pi]

# fixed pressure boundary
    

# height 
# h(t,x) = h(x)

# h(x) = h0 + delta * sin(k0 x) * cos(k1 z)
h0, h_delta = [0.2, 0.15] # delta < h0 so height positive
h_alpha, h_beta = [3, 2]
str_h = "%0.2f + %0.2f \sin(%d x) \cos(%d z)"%(h0, h_delta, h_alpha, h_beta) #for graph title

def h(t, X):
    return h0 + h_delta * np.sin(h_alpha*X[0]) * np.cos(h_beta*X[1])

def h_dX(t, X): # = [h_x, h_z]
    x, z = X
    return [h_delta * h_alpha * np.cos(h_alpha*x)* np.cos(h_beta*z), -h_delta * h_beta * np.sin(h_alpha*x) * np.sin(h_beta*z)]

def eval_height(ts, xs, zs): 
    hs = np.zeros((len(ts), len(xs), len(zs)))
    for k in range(len(ts)): 
        for i in range(len(xs)):
            for j in range(len(zs)):
                hs[k][i][j] = h(ts[k], (xs[i], zs[j]))
    return hs
    

#------------------------------------------------------------------------------
# Manf. solution & RHS
#------------------------------------------------------------------------------
# (h^3 P')' = f(h)

# Pressure P(xy)

p_alpha, p_beta = [2,-1]
def p(X):
    return -np.sin(p_alpha*X[0])*np.cos(p_beta*X[1])

str_p = "-\sin(%d x)\cos(%d z)"%(p_alpha, p_beta) #for graph title

def eval_p(Xs):
    p_exact = np.zeros(len(Xs))
    for i in range(len(Xs)):
        p_exact[i] = p(Xs[i])
    return p_exact

def p_dX(X): 
    p_x = -p_alpha * np.cos(p_alpha*X[0]) * np.cos(p_beta*X[1])
    p_z = p_beta * np.sin(p_alpha*X[0]) * np.sin(p_beta*X[1])
    return [p_x, p_z]

def p_dXX(X):
    p_xx = p_alpha**2 * np.sin(p_alpha*X[0]) * np.cos(p_beta*X[1])
    p_zz = p_beta**2 * np.sin(p_alpha*X[0]) * np.cos(p_beta*X[1])
    return [p_xx , p_zz]

# calculate RHS f(x,t)

def f(t,X): # = h^3 p'' + (h^3)' p' 
    p_x, p_z = p_dX(X)
    p_xx, p_zz = p_dXX(X)
    h_x, h_z = h_dX(t, X)
    
    return h(t,X)**3 * (p_xx + p_zz) + 3*h(t,X)**2 * (h_x*p_x + h_z*p_z)

def eval_f(ts, Xs):
    f_n = np.zeros((len(ts), len(Xs))) #f[t][xz]
    for i in range (len(ts)):
        for j in range(len(Xs)):
            f_n[i][j] = f(ts[i], Xs[j])
    return f_n
    
#------------------------------------------------------------------------------
# numerical solution
#------------------------------------------------------------------------------   
# plot settings
view_theta = 30 # pan view angle up/down
view_phi = 45  # pan view angle left-right
# figs = [0: return D p_n, 1: slice (x, z0), 2: slice (x0, z), 3: plot 3D]
# BCs = [0: periodic, 1: fixed pressure]

def solve(Nt, Nx, Nz, fig=3, BC=0):
    start = perf_counter()
    
    # grid 
    dt = (tf - t0)/Nt
    dx = (xb - xa)/Nx
    dz = (zb - za)/Nz

    ts = [t0 + i*dt for i in range(Nt)]
    xs = [xa + i*dx for i in range(Nx)]
    zs = [za + i*dz for i in range(Nz)]
    
    # flat space [(x0,z0), (x0, z1), ..., (x0, zN), (x1, z0), ..., (xN, zN)]
    Xs = np.zeros((Nx*Nz,2)) 
    for i in range(len(Xs)):
        Xs[i] = [xs[i//Nz], zs[i%Nz]]  

    # height h(t,X) = h[t][x][z]
    hs = eval_height(ts, xs, zs)
    
    # RHS f(t,X) = f_n[t][xz]
    f_n = eval_f(ts, Xs)
        
    # exact solution p(t, X) = p[xzs]
    p_exact = eval_p(Xs)
    
    # Loop over time t = ts[i]
    inf_norms = np.zeros(Nt) 
    for k in range(Nt): #time: ts[k]
    
        # construct finite difference matrix D
        D = np.zeros((Nx*Nz, Nx*Nz))
        
        for i in range(Nx): # xs[i]
            
            # Build rows [i*Nz : (i+1)*Nz]
            
            # Initialize 5 diagonals... 
            P_center = np.zeros(Nz)  #P(i,j)    central diagonal 
            P_west = np.zeros(Nz)    #P(i, j-1) lower diagonal
            P_east = np.ones(Nz)     #P(i, j+1) upper diagonal 
            P_south = np.ones(Nz)    #P(i-1, j) left diagonal
            P_north = np.ones(Nz)    #P(i+1, j) north diagonal
        
            for j in range(Nz): # ys[j]
            
                # Find h(t,X) on grid                  
                h_c = hs[k][i][j]        #h(i,j)    
                h_N = hs[k][(i+1)%Nx][j] #h(i+1, j)
                h_S = hs[k][(i-1)%Nx][j] #h(i-1, j)
                h_E = hs[k][i][(j+1)%Nz] #h(i, j+1)
                h_W = hs[k][i][(j-1)%Nz] #h(i, j-1)
                
                #P(i,j) central diagonal
                P_center_x = -(h_N**3 + 3*h_c*h_N**2 + 3*h_N*h_c**2 + 2*h_c**3 + 3*h_S*h_c**2 + 3*h_c*h_S**2 + h_S**3)/(8* dx**2) #Pij from Px
                P_center_z = -(h_E**3 + 3*h_c*h_E**2 + 3*h_E*h_c**2 + 2*h_c**3 + 3*h_W*h_c**2 + 3*h_c*h_W**2 + h_W**3)/(8* dz**2) #Pij from Py
                P_center[j] = P_center_x + P_center_z
                
                #P(i+1, j) right diagonal
                P_north[j] = (h_c**3 + 3*h_c*h_N**2 + 3*h_N*h_c**2 + h_N**3)/(8*dx**2)
                
                #P(i-1, j) left diagonal
                P_south[j] = (h_c**3 + 3*h_c*h_S**2 + 3*h_S*h_c**2 + h_S**3)/(8*dx**2)
                
                #P(i, j-1) lower diagonal
                P_west[j]  = (h_c**3 + 3*h_c*h_W**2 + 3*h_W*h_c**2 + h_W**3)/(8*dz**2)
                
                #P(i, j+1) upper diagonal 
                P_east[j]  = (h_c**3 + 3*h_c*h_E**2 + 3*h_E*h_c**2 + h_E**3)/(8*dz**2)


            if BC == 0:
                
                # Make three (Ny * Ny) matrices with these diagonals
                D_i_left = np.diagflat(P_south)
                D_i_center = np.diagflat(P_west[1:Nz], -1) + np.diagflat(P_center) + np.diagflat(P_east[0:Nz-1], 1)
                D_i_right = np.diagflat(P_north)
                
                # add corners to D_i_center
                D_i_center[0][Nz-1] = P_west[0]
                D_i_center[Nz-1][0] = P_east[Nz-1]
                
                # input rows [i*Nx : (i+1)*Nx] into D
                if i == 0: # top: [ center | right | ... | left ]
                    D[0:Nz, 0:Nz] = D_i_center
                    D[0:Nz, Nz:2*Nz] = D_i_right
                    D[0:Nz, (Nx-1)*Nz:Nx*Nz] = D_i_left
                elif i == Nx-1: #bottom: [ right | ... | left | center ]
                    D[(Nx-1)*Nz:Nx*Nz, 0:Nz] = D_i_right
                    D[(Nx-1)*Nz:Nx*Nz, (Nx-2)*Nz:(Nx-1)*Nz] = D_i_left
                    D[(Nx-1)*Nz:Nx*Nz, (Nx-1)*Nz:Nx*Nz] = D_i_center
                else: #standard: [ ... | left | center | right | ... ]
                    D[i*Nz:(i+1)*Nz, (i-1)*Nz:i*Nz] = D_i_left
                    D[i*Nz:(i+1)*Nz, i*Nz:(i+1)*Nz] = D_i_center
                    D[i*Nz:(i+1)*Nz, (i+1)*Nz:(i+2)*Nz] = D_i_right
            
            elif BC == 1:
                p0 = 0
                
                if i == 0:
                    P_center = np.ones(Nz)
                    P_west = np.ones(Nz)
                    P_east = np.ones(Nz)
                    
                    P_north[0] = 1
                    P_north[Nz-1] = 1
                    
                    P_south = np.ones(Nz)
                    
                    f_n[k][0:Nz] = p0 #p(xa, zs[0:Nz])
                    
                elif i == 1:
                    P_south = np.ones(Nz)
                    
                    P_north[0] = 1
                    P_north[Nz-1] = 1
                    
                    P_center[0] = 1
                    P_center[Nz-1] = 1
                    P_west[0] = 1
                    P_west[1] = 1
                    P_east[Nz-2] = 1
                    P_east[Nz-1] = 1
                    
                    f_n[k][Nz] = p0 #p(xs[1], za)
                    f_n[k][2*Nz - 1] = p0 #p(xs[1], zb)
                    

                elif i == Nx-2:
                    P_south[0] = 1
                    P_south[Nz-1] = 1
                    
                    P_center[0] = 1
                    P_center[Nz-1] = 1
                    P_west[0] = 1
                    P_west[1] = 1
                    P_east[Nz-2] = 1
                    P_east[Nz-1] = 1
                    
                    P_north = np.ones(Nz)
                    
                    f_n[k][(Nx-2)*Nz] = p0# p(xs[Nx-2], za)
                    f_n[k][Nx*(Nz-1)-1] = p0# p(xs[Nx-2], zb)
                    
                
                elif i == Nx-1:
                    P_north = np.ones(Nz)
                    
                    P_south[0] = 1
                    P_south[Nz-1] = 1
                    
                    P_center = np.ones(Nz)
                    P_west = np.ones(Nz)
                    P_east = np.ones(Nz)
                    
                    f_n[k][(Nx-1)*Nz] = p0# p(xb, za)
                    f_n[k][(Nx-1)*Nz+1:Nx*Nz - 2] = p0#p(xb, zs[1:Nx-1])
                    f_n[k][Nx*Nz - 1] = p0#p(xb, zb)
                    
        
                else: 
                    
                    P_south[0] = 1
                    P_south[Nz-1] = 1
                    
                    P_north[0] = 1
                    P_north[Nz-1] = 1
                    
                    P_center[0] = 1
                    P_center[Nz-1] = 1
                    P_west[0] = 1
                    P_west[1] = 1
                    P_east[Nz-2] = 1
                    P_east[Nz-1] = 1
                    
                    f_n[k][i*Nz] = p0#p(xs[i], zs[0]) 
                    f_n[k][Nx * (i+1) - 1] = p0#p(xs[Nx * (i+1) - 1], zs[Nz-1])
                
                # Make three (Ny * Ny) matrices with these diagonals
                D_i_left = np.diagflat(P_south)
                D_i_center = np.diagflat(P_west[1:Nz], -1) + np.diagflat(P_center) + np.diagflat(P_east[0:Nz-1], 1)
                D_i_right = np.diagflat(P_north)
                
                # add corners to D_i_center
                D_i_center[0][Nz-1] = P_west[0]
                D_i_center[Nz-1][0] = P_east[Nz-1]
                
                # input rows [i*Nx : (i+1)*Nx] into D
                if i == 0: # top: [ center | right | ... | left ]
                    D[0:Nz, 0:Nz] = D_i_center
                    D[0:Nz, Nz:2*Nz] = D_i_right
                    D[0:Nz, (Nx-1)*Nz:Nx*Nz] = D_i_left
                elif i == Nx-1: #bottom: [ right | ... | left | center ]
                    D[(Nx-1)*Nz:Nx*Nz, 0:Nz] = D_i_right
                    D[(Nx-1)*Nz:Nx*Nz, (Nx-2)*Nz:(Nx-1)*Nz] = D_i_left
                    D[(Nx-1)*Nz:Nx*Nz, (Nx-1)*Nz:Nx*Nz] = D_i_center
                else: #standard: [ ... | left | center | right | ... ]
                    D[i*Nz:(i+1)*Nz, (i-1)*Nz:i*Nz] = D_i_left
                    D[i*Nz:(i+1)*Nz, i*Nz:(i+1)*Nz] = D_i_center
                    D[i*Nz:(i+1)*Nz, (i+1)*Nz:(i+2)*Nz] = D_i_right
                    
                
                

        # closure: assume sum p'' = 0 
        D[Nx*Nz - 1] = 1 # [1 ... 1]
        f_n[k][Nx*Nz - 1] = 0
        
        # solve for p_n = [p00, p01, p02, ..., p0N, p10, p11, ..., pNN]
        solve_start = perf_counter()
        p_n = spl.spsolve(D, f_n[k])
        
        #p_n = np.linalg.solve(D, f_n[k])
        solve_end = perf_counter()
        
        
        # Plotting at t=ts[0]
        if k == 0: 
        
            if fig == 1: #(x, z0, p)
                z0 = za
                p_n_1D_x = np.zeros(Nx)
                p_exact_1D_x = np.zeros(Nx)

                for i in range(Nx):
                    p_n_1D_x[i] = p_n[i*Nz + z0]
                    p_exact_1D_x[i] = p_exact[i*Nz+z0]
                    
                pp.figure()
                pp.plot(xs, p_n_1D_x, label="$p_n(y_0)$")
                pp.plot(xs, p_exact_1D_x, label="$p_n(y_0)$")
                pp.title("$p= %s$ | $h=%s$ | $N_x=%d$"%(str_p, str_h, Nx))
                pp.xlabel('x')
                pp.ylabel('p')
                pp.legend()

            elif fig == 2: #(x0, z, p)
                x0 = xa
                p_n_1D_z = p_n[x0*Nz:x0*Nz + Nz]
                p_exact_1D_z = p_exact[x0*Nz:x0*Nz + Nz]
                    
                pp.figure()
                pp.plot(zs, p_n_1D_z, label="$p_n(x_0)$")
                pp.plot(zs, p_exact_1D_z, label="$p(x_0)$")
                pp.title("$p= %s$ | $h=%s$ | $N_y=%d$"%(str_p, str_h, Nz))
                pp.xlabel('z')
                pp.ylabel('p')
                pp.legend()
        
            elif fig == 3: #(x, z, p)
                p_n_2D = np.zeros((Nx, Nz))
                p_exact_2D = np.zeros((Nx, Nz))
                
                
                for i in range(Nx):
                    for j in range(Nz):
                        
                        p_n_2D[i,j] = p_n[i*Nz + j]
                        p_exact_2D[i,j] = p_exact[i*Nz + j]
                        
                X, Y = np.meshgrid(xs, zs)        
                    
                pp.figure()
                ax = pp.axes(projection='3d')
                ax.plot_surface(X.T, Y.T, p_n_2D, label="$p_N$", rstride=1, cstride=1,cmap='viridis')
                pp.title("$p= %s$ | $h=%s$ | $N_x=%d$"%(str_p, str_h, Nx))
                pp.xlabel('x')
                pp.ylabel('y')
                ax.view_init(view_theta, view_phi)
                
                pp.figure()
                ax = pp.axes(projection='3d')
                ax.plot_surface(X.T, Y.T, p_exact_2D, label="$p_exact$", rstride=1, cstride=1, cmap='viridis')
                pp.title("$p = %s$ | $h=%s$ | exact"%(str_p, str_h))
                pp.xlabel('x')
                pp.ylabel('z')
                ax.view_init(view_theta, view_phi)
    
        # error at t = ts[k]     
        inf_norms[k] = np.max(np.abs(np.subtract(p_exact,p_n)))
        
        print("Solved Nx=%d, Nz=%d"%(Nx, Nz))    
        print("Error %0.5f at t=%0.2f"%(inf_norms[k], ts[0]))
        
        end = perf_counter()
        
        print("Total run time = %0.3f"%(end-start))
        print("Solve time = %0.3f"%(solve_end-solve_start))
        #print("solve/run ratio %0.4f"%((solve_end-solve_start)/(end-start)))
        print("\n")
    return inf_norms

#------------------------------------------------------------------------------
# Convergence
#------------------------------------------------------------------------------   
def conveg(trials=6, N0=5):
    Nx = N0
    Nz = N0
    Nt = 1
    
    infNorms = np.zeros(trials)
    dxs = np.zeros(trials)
    dxs_sqr = np.zeros(trials)
    
    for i in range(trials):
        
        infNorms[i] = np.max(solve(Nx, Nz, Nt, 0)) # max over time
        dxs[i] = ((xb - xa) / Nx)
        dxs_sqr[i] = dxs[i]**2
        
        Nx += 50
        Nz += 50
    
    pp.figure(trials+1)
    pp.loglog(dxs, dxs_sqr, color='r', label='$dx^2$')
    pp.loglog(dxs, infNorms, color='b', label='Linf Error')
    
    pp.xlabel('dx')

    pp.legend()
    pp.title('$N_x$=%i to $N_x$=%i with %i trials'%(N0, Nx-50, trials))

