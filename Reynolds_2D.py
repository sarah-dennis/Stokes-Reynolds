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
import _graphics

#------------------------------------------------------------------------------
# Domain
#------------------------------------------------------------------------------

# -- time --------------------------------------------------------------------
t0, tf = [0, 0]

# -- space --------------------------------------------------------------------
xa, xb = [0, 6]

za, zb = [0, 6]

  
# flat space [(x0,z0), (x0, z1), ..., (x0, zN), (x1, z0), ..., (xN, zN)]
def grid_eval_time(f, ts, xs, zs):
    f_n = np.zeros((len(ts), len(xs)*len(zs))) #f[t][xz]
    for k in range (len(ts)):
        for i in range(len(xs)):
            for j in range(len(zs)):
                f_n[k][i*len(zs)+j] = f(ts[k], xs[i], zs[j])
    return f_n

def grid_eval_notime(f, xs, zs):
    f_n = np.zeros(len(xs)*len(zs)) #f[t][xz]
    for i in range(len(xs)):
        for j in range(len(zs)):
            f_n[i*len(zs)+j] = f(xs[i], zs[j])
    return f_n

# -- height --------------------------------------------------------------------

# option 1: sinusoidal height
## h(X) = h0 + delta * sin(k0 x) * cos(k1 z)

# h0, h_delta = [0.2, 0.15] # delta < h0 
# h_alpha, h_beta = [1, 1]

# str_h = "%0.2f + %0.2f \sin(%d x) \cos(%d z)"%(h0, h_delta, h_alpha, h_beta) #for graph title

# def h(t, x, z):
#     return h0 + h_delta * np.sin(h_alpha*x) * np.cos(h_beta*z)

# def h_dX(t, x, z): # = [h_x, h_z]
#     return [h_delta * h_alpha * np.cos(h_alpha*x)* np.cos(h_beta*z), -h_delta * h_beta * np.sin(h_alpha*x) * np.sin(h_beta*z)]

# option 2: wedge height (must use BC = 1 or 2)
## h(X) = hf + mx

h0, hf = [0.1, 0.01] # h0 > hf > 0
h_m = (2*hf - h0)/(xb - xa)

str_h = "%0.2f + %0.2f x"%(h0, h_m) #for graph title

def h(t, x, z):
    return hf + h_m*x

def h_dX(t, x, z): # = [h_x, h_z]
    return [h_m, 0]

# h(t) vs x(t) -- > which is really h(x(t))
def h_dt(t, x, z):
    return 0


#------------------------------------------------------------------------------
# Manf. solution
#------------------------------------------------------------------------------
# (h^3 P')' = f(h)

# Pressure P(xy)

p_alpha, p_beta = [1,1]
def p(x, z):
    return -np.sin(p_alpha*x)*np.cos(p_beta*z)

str_p = "p(x, z) = -\sin(%d x)\cos(%d z)"%(p_alpha, p_beta) #for graph title


def p_dX(x, z): 
    p_x = -p_alpha * np.cos(p_alpha*x) * np.cos(p_beta*z)
    p_z = p_beta * np.sin(p_alpha*x) * np.sin(p_beta*z)
    return [p_x, p_z]

def p_dXX(x, z):
    p_xx = p_alpha**2 * np.sin(p_alpha*x) * np.cos(p_beta*z)
    p_zz = p_beta**2 * np.sin(p_alpha*x) * np.cos(p_beta*z)
    return [p_xx , p_zz]

def plot_pressure():
    Nx = _graphics.fig_Nx
    Nz = _graphics.fig_Nz
    dx = (xb - xa) / (Nx-1)
    dz = (zb - za) / (Nz-1)
    xs = [xa + k*dx for k in range(Nx)]
    zs = [za + k*dz for k in range(Nz)]
    p_2D = np.zeros((Nx, Nz))
    for i in range(Nx):
        for j in range(Nz):
            p_2D[i][j] = p(xs[i], zs[j])
    _graphics.plot_3D(p_2D, xs, zs, str_p)
    
#------------------------------------------------------------------------------
# RHS f(h(x,t)) = ht + 0.5 ((h Ux)x + (h Uz)z)
#------------------------------------------------------------------------------

# calculate RHS f(x,t) = ht + 0.5 (hx Ux + h Uxx + hz Uz + h Uzz)

# manufactured RHS
# f = h^3 p'' + (h^3)' p' 
def f(t, x, z): 
    p_x, p_z = p_dX(x, z)
    p_xx, p_zz = p_dXX(x, z)
    h_x, h_z = h_dX(t, x, z)
    
    return h(t,x, z)**3 * (p_xx + p_zz) + 3*h(t,x, z)**2 * (h_x*p_x + h_z*p_z)

# Reynolds RHS ?
# def f(t, X):
#     h_x, h_z = h_dX(t, X)
#     h_t = h_dt(t,X)
#     h_ = h(t,X)

#     return  0.5*h_*(h_x + h_z)


#------------------------------------------------------------------------------
# numerical solution
#------------------------------------------------------------------------------   
# figs = [0: none, 1: numerical pressure, 2: exact pressure, 3: error]

# BCs = [0: periodic, 1: fixed pressure, 2: fixed x & periodic z]
BCs_txt = ['x & z periodic', 'x & z prescribed', 'x prescribed, z periodic']

def solve(Nt=1, Nx=100, Nz=100, fig=3, BC=1):
    start = perf_counter()
    
    # grid 
    dt = 0
    dx = (xb - xa)/(Nx-1)
    dz = (zb - za)/(Nz-1)

    ts = [t0 + i*dt for i in range(Nt)]
    xs = [xa + i*dx for i in range(Nx)]
    zs = [za + i*dz for i in range(Nz)]
  
    # height h(t,X) = h[t][xz]
    # hs = grid_eval_time(f, ts, xs, zs)
    
    # RHS f(t,X) = f_n[t][xz]
    f_n = grid_eval_time(f, ts, xs, zs)
        
    # exact solution p(t, x, z) = p[xzs]
    p_exact = grid_eval_notime(p, xs, zs)
    
    # Loop over time t = ts[i]
    errs = np.zeros(Nt) 
    
    for k in range(Nt): #time: ts[k]
        t = ts[k]
        # construct finite difference matrix D
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
                h_c = h(t, xs[i], zs[j])         #h(i,j)    
                h_N = h(t, xs[(i+1)%Nx], zs[j])  #h(i+1, j)
                h_S = h(t, xs[(i-1)%Nx],  zs[j]) #h(i-1, j)
                h_E = h(t, xs[i], zs[(j+1)%Nz])  #h(i, j+1)
                h_W = h(t, xs[i], zs[(j-1)%Nz])  #h(i, j-1)

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
                f_n[k][Nx*Nz - 1] = 0
            
            elif BC == 1: #fixed in x and z
            
                #prescribe pressure at x0 and xf
                if i == 0 or i == Nx-1: 
                
                    P_center = np.ones(Nz)
                    
                    P_west = np.zeros(Nz)   
                    P_east = np.zeros(Nz)     
                    P_south = np.zeros(Nz)    
                    P_north = np.zeros(Nz)
                    
                    f_n[k][i*Nz:(i+1)*Nz] = [p(xs[i], zs[j]) for j in range(0, Nz)] #p(xs[i], zs[0:Nz])
            
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
                    
                    f_n[k][i*Nz] = p(xs[i], zs[0])
                    f_n[k][i*Nz + Nz-1] = p(xs[i], zs[Nz-1])

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
                    
                    f_n[k][i*Nz:(i+1)*Nz] = [p(xs[i], zs[j]) for j in range(0, Nz)] #p(xa, zs[0:Nz])
               
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
        p_n = spl.spsolve(D, f_n[k])
        #p_n = np.linalg.solve(D, f_n[k])
        solve_end = perf_counter()
        
        # Plotting at t0 = ts[k=0]
        if k == 0: 
            if fig == 1: # p_n(x, z)
                p_n_2D = np.zeros((Nx, Nz))
                p_exact_2D = np.zeros((Nx, Nz))
                
                for i in range(Nx):
                    for j in range(Nz):
                        p_n_2D[i,j] = p_n[i*Nz + j]
                        p_exact_2D[i,j] = p_exact[i*Nz + j]
                        
                title = "$p= %s$ | $h=%s$ \n $N_x, N_z=%d, %d$ | BC: %s"%(str_p, str_h, Nx, Nz, BCs_txt[BC])
                _graphics.plot_3D(p_n_2D, xs, zs, title)
            
            elif fig==2: # p_exact(x, z)
                p_exact_2D = np.zeros((Nx, Nz))
                
                for i in range(Nx):
                    for j in range(Nz):
                        p_exact_2D[i,j] = p_exact[i*Nz + j]
                        
                title = "$p = %s$ | $h=%s$ | exact"%(str_p, str_h)
                _graphics.plot_3D(p_exact_2D, xs, zs, title)
            
            elif fig==3: # Error
            
                p_n_2D = np.zeros((Nx, Nz))
                p_exact_2D = np.zeros((Nx, Nz))
                
                for i in range(Nx):
                    for j in range(Nz):
                        p_n_2D[i,j] = p_n[i*Nz + j]
                        p_exact_2D[i,j] = p_exact[i*Nz + j]
                
                title = "Error | $N_x=%d$ , $N_z=%d$ | $dx =%.3f$ $dz = %.3f$ \n | $h(x,z) = %s$ | BC: %s"%(Nx, Nz, dx, dz, str_h, BCs_txt[BC])
                _graphics.plot_3D(p_n_2D-p_exact_2D, xs, zs, title)
               
    
        # error at t = ts[k]     
        errs[k] = np.max(np.abs(np.subtract(p_exact,p_n)))
        
        print("Solved x:[%.1f,%.1f] z:[%.1f, %.1f] with (Nx=%d, Nz=%d)~(dx=%.3f, dz=%.3f) and BC: %s"%(xa, xb, za, zb, Nx, Nz, dx, dz, BCs_txt[BC]))    
        print("Error %0.5f at t=%0.2f"%(errs[k], ts[0]))
        
        end = perf_counter()
        
        print("Solve time = %0.3f"%(solve_end-solve_start))
        print("Run time = %0.3f"%(end-start))
        #print("solve/run ratio %0.4f"%((solve_end-solve_start)/(end-start)))
        print("\n")
    return errs

#------------------------------------------------------------------------------
# Convergence -- fix 
#------------------------------------------------------------------------------ 
# BCs = [0: periodic, 1: fixed pressure, 2: fixed x & periodic z]
def conveg(trials=5, N0=10, BC=1):
    Nx = N0
    Nz = N0
    Nt = 1
    
    dNx = 2
    dNz = 2
    
    infNorms = np.zeros(trials)
    
    dxs = np.zeros(trials)
    dxs_sqr = np.zeros(trials)
    
    fig = 0

    
    for i in range(trials):
        
        if i == trials-1: #make 3D plot at the last iteration
            fig = 4
       
        print('starting trial %d of %d'%(i+1, trials))
            
        infNorms[i] = np.max(solve(Nt, Nx, Nx, fig, BC)) #max over time
        
        dx = (xb - xa) / (Nx-1)
        dz = (zb - za) / (Nz-1)
        dxs[i] = dx
        dxs_sqr[i] = (dx**2) + (dz**2)
        
        Nx *= dNx
        Nz *= dNz
    
    pp.figure(trials+1)
    pp.loglog(dxs, dxs_sqr, color='r', label='$dx^2 + dz^2$')
    pp.loglog(dxs, infNorms, color='b', label='$L_\inf$ Error')
    pp.xticks(dxs[::2])
    pp.xlabel('$dx$')

    pp.legend()
    pp.title('$N_x$=%i to $N_x$=%i over %i trials | Bc: %s \n $h(x,z)=%s$ | $p(x,z)=%s$'%(N0, Nx-dNx, trials, BCs_txt[BC], str_h, str_p))
