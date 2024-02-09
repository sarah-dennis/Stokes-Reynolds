#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 22 14:59:46 2024

@author: sarahdennis
"""
import numpy as np
import csv

from scipy.sparse.linalg import LinearOperator
from scipy.sparse.linalg import bicgstab

import graphics as graph



x0 = 0
xL = 1

y0 = 0
yL = 2

slope = yL/(xL/2)
U = 1

Re = 1

N = 500
h = 1/N

n = xL*N + 1
m = yL*N + 1

stab_tol = 1e-5
trials = 10
plot_every = 1

name = "guptaKalita_N%d.csv"%N 

def run_new(filename):
    u = np.ones(n*m)
    v = np.ones(n*m)
    
    u, v = run(u,v)
    
    with open(filename, 'w', newline='') as file:
        writer = csv.writer(file, delimiter=' ')
        for i in range(n*m):
            writer.writerow([u[i], v[i]])
    
def run_load(filename):

    with open(filename, newline='') as file:
        reader = csv.reader(file)
        u = np.ones(n*m)
        v = np.ones(n*m)
        
        for i in range(n*m):
            ui, vi = next(reader)[0].split(' ')
            u[i] = float(ui)
            v[i] = float(vi)

    u, v = run(u,v)
    
    with open(filename, 'w', newline='') as file:
        writer = csv.writer(file, delimiter=' ')
        for i in range(n*m):
            writer.writerow([u[i], v[i]])
    
    
def run(u,v):
    rhs = make_rhs(u,v)
    M = Dpsi_linOp()

    for i in range(trials): 
        print("trial %d of %d"%(i+1, trials))
        rhs = make_rhs(u, v)
        stream, exit_flag = bicgstab(M, rhs, atol=stab_tol)
        u, v = uv_approx(u, v, stream)
        
        if i % plot_every == 0:
            make_plots(u, v, i)
            
    return u, v
    
# -----------------------------------------------------------------------------
# rhs = c0*A - c1*(B - C)    

# A = u_S - u_N + v_E - v_W
# B = v_C * (u_E + u_W + u_N + u_S)
# C = u_C * (v_E + v_W + v_N + v_S)

c0 = 3*h
c1 = 0.5 * h**2 * Re
def make_rhs(u, v): #
    rhs = np.zeros(m*n)
    for k in range(n*m):
        i = k % n
        j = k//n
        
        x = h*i
        y = h*j
       
        # boundary
        if i == 0 or j == 0 or i == n-1 or j == m-1:
            rhs[k] = 0
          
        # dead zones 
        elif y <= -slope * x + yL or y <= slope * x - yL:
            rhs[k] = 0
        
        # interior
        else:
            # (u,v) at 9 point stencil
            #k = j*n + i
            u_C = u[k]
            v_C = v[k]
            
            # North (i, j+1)
            if j+1 == m-1: 
                u_N = U
                v_N = 0
            else:
                u_N = u[(j+1)*n + i]
                v_N = v[(j+1)*n + i]
                
            # East (i+1, j)
            x_E = h*(i+1)
            if i+1 == n-1 or y <= slope * x_E - yL: 
                u_E = 0
                v_E = 0
            else:
                v_E = v[j*n + i+1]
                u_E = u[j*n + i+1]
            
            # South (i, j-1)
            y_S = h*(j-1)
            if j-1 == 0 or y_S <= slope*x - yL or y_S <= -slope*x + yL: 
                v_S = 0
                u_S = 0
            else:
                v_S = v[(j-1)*n + i]
                u_S = u[(j-1)*n + i]
                
            # West (i-1, j)  
            x_W = (i-1)*h
            if i-1 == 0 or y <= -slope * x_W - yL: 
                v_W = 0
                u_W = 0
            else:
                v_W = v[j*n + (i-1)]
                u_W = v[j*n + (i-1)]
                
            A = u_S - u_N + v_E - v_W
            B = v_C * (u_E + u_W + u_N + u_S)
            C = u_C * (v_E + v_W + v_N + v_S)

            rhs[k] = c0 * A + c1 * (B - C)
            
    return rhs

# -----------------------------------------------------------------------------

# new_u[k] = c2 * (psi_N - psi_S) - c3 * (u_N + u_S)
# new_v[k] = -c2 * (psi_E - psi_W) - c3 * (v_E + v_W)

c2 = 3/(4*h)
c3 = 1/4
def uv_approx(u, v, psi):
    new_u = np.zeros(n*m)
    new_v = np.zeros(n*m)

    for k in range(n*m):
        i = k % n
        j = k//n
        
        x = h*i
        y = h*j
    
        # y=yL boundary
        if j == m-1: 
            new_u[k] = U
            new_v[k] = 0 
            
        # other side boundaries
        elif i == 0 or j == 0 or i == n-1:
            new_u[k] = 0
            new_v[k] = 0 
            
        # dead zones
        elif y <= -slope * x + yL or y <= slope * x - yL:
            
            new_u[k] = 0
            new_v[k] = 0 
                
        else: #interior
            
            # (u,v) at 4 point stencil

            # North (i, j+1)
            if j+1 == m-1: 
                u_N = U
                psi_N = 0
            else:
                u_N = u[(j+1)*n + i]
                psi_N = psi[(j+1)*n + i]
                
            # East (i+1, j)
            x_E = h*(i+1)
            if i+1 == n-1 or y <= slope * x_E - yL: 
                v_E = 0
                psi_E = 0
            else:
                v_E = v[j*n + i+1]
                psi_E = psi[j*n + i]
            
            # South (i, j-1)
            y_S = h*(j-1)
            if j-1 == 0 or y_S <= slope*x - yL or y_S <= -slope*x + yL: 
                u_S = 0
                psi_S = 0
            else:
                u_S = u[(j-1)*n + i]
                psi_S = psi[(j-1)*n + i]
                
            # West (i-1, j)  
            x_W = h*(i-1)
            if i-1 == 0 or y <= -slope * x_W - yL: 
                v_W = 0
                psi_W = 0
            else:
                v_W = v[j*n + (i-1)]
                psi_W = psi[j*n + i-1]
   
            new_u[k] = c2 * (psi_N - psi_S) - c3 * (u_N + u_S)
            new_v[k] = -c2 * (psi_E - psi_W) - c3 * (v_E + v_W)
    
    return new_u, new_v

# -----------------------------------------------------------------------------

class Dpsi_linOp(LinearOperator):

    def __init__(self):
        #  globals: Nx, Ny, h, slope, yL
        self.shape = ((n*m, n*m))        
        self.dtype = np.dtype('f8')
        self.mv = np.zeros(n*m) 
        
    def _matvec(self, v): # v:= psi[i*n + j]  

        for k in range(n*m):
            i = k%n
            j = k//n
            
            x = h*i
            y = h*j

            # outer boundary
            if i == 0 or i == n-1 or j == 0 or j == m-1:
                self.mv[k] = 0
    
            # dead zones: 
            elif y <= -slope * x + yL or y <= slope * x - yL:
                self.mv[k] = 0

            #interior   
            else: 
                # psi[k] at 9 point stencil         
                
                v_C = v[k]
        
                #North (i, j+1)
                y_N = h*(j+1)
                if j+1 == m-1: 
                    v_N = 0
                else:
                    v_N = v[(j+1)*n + i]
                    
                #East (i+1, j)
                x_E = h*(i+1)
                if i+1 == n-1 or y <= slope * x_E - yL:
                    v_E = 0
                else:
                    v_E = v[j*n + (i+1)]
                        
                #South (i, j-1)
                y_S = h*(j-1)
                if j-1 == 0 or y_S <= slope*x - yL or y_S <= -slope*x + yL: 
                    v_S = 0
                else:
                    v_S = v[(j-1)*n + i]
                
                #West (i-1,j)
                x_W = h*(i-1)
                if i-1 == 0 or y <= -slope * x_W - yL: 
                    v_W = 0
                else:
                    v_W = v[j*n + (i-1)]
                    
                #NorthEast (i+1, j+1)
                if i+1 == n-1 or j+1 == m-1 or y_N <= slope * x_E - yL: 
                    v_NE = 0
                else:
                    v_NE  = v[(j+1)*n + (i+1)]
                
                #NorthWest (i-1, j+1)
                if i-1 == 0 or j+1 == m-1 or y_N <= -slope * x_W - yL:
                    v_NW = 0
                else:
                    v_NW = v[(j+1)*n + (i-1)]
                
                #SouthEast (i+1, j-1)
                if i+1 == n-1 or j-1 == 0 or y_S <= slope*x_E - yL or y_S <= -slope*x_E + yL: 
                    v_SE = 0
                else:
                    v_SE = v[(j-1)*n + i+1]
                
                #SouthWest (i-1, j-1)
                if i-1 == 0 or j-1 == 0 or y_S <= slope*x_W - yL or y_S <= -slope*x_W + yL:
                    v_SW = 0
                else:
                    v_SW = v[(j-1)*n + i-1]

        
                self.mv[k] = v_SW - 8*v_S + v_SE + -8*v_W + 28*v_C -8*v_E + v_NW -8*v_N + v_NE

        return self.mv

#-----------------------------------------------------------------------------

def make_plots(u, v, k):
# Grid domain
    xs = np.linspace(x0, xL, n)
    ys = np.linspace(y0, yL, m)

# velocity 
    u_2D = u.reshape((m,n))
    v_2D = v.reshape((m,n))
    ax_labels = ['x', 'y']
    title = 'velocity stream-plot ($N=%d$, $trials=%d/%d$)'%(N, k, trials)
    graph.plot_stream(u_2D, v_2D, xs, ys, title, ax_labels)

# vorticity 

    uy_2D = np.gradient(u_2D, h, axis=0)
    vx_2D = np.gradient(v_2D, h, axis=1)
    w = np.zeros((m,n))
    for j in range(m):
        for i in range(n):   
            w[j,i] = vx_2D[j,i] - uy_2D[j,i]
            if w[j,i] == 0:
                w[j,i] = None

    ax_labels = ['$\omega = v_x - u_y$', 'x', 'y']
    title = 'vorticity stream-plot ($N=%d$, $trials=%d/%d$)'%(N, k, trials)
    graph.plot_heatMap(w, xs, ys, title, ax_labels)


