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

stab_tol = 1e-5
n_trials= 3
plot_mod = 1

class triangle():
    def __init__(self, x0, xL, y0, yL, U, Re, N):
        self.x0 = x0
        self.xL = xL
        self.y0 = y0
        self.yL = yL
        self.slope = yL/(xL/2)
        self.U = U
        self.Re = Re
        self.N = N
        self.h = 1/N
        self.Nx = xL*N + 1
        self.Ny = yL*N + 1
        self.filename = "guptaKalita_N%d.csv"%N

class biswasKalita(triangle):
    def __init__(self, N):
        x0 = 0
        xL = 1
        y0 = 0
        yL = 2

        U = 1
        Re = 1
        super().__init__(x0, xL, y0, yL, U, Re, N)


#Run with U,V = 1 and write to csv
def run_new(tri):
    n = tri.Nx
    m = tri.Ny
    
    u = np.ones(n*m)
    v = np.ones(n*m)
    
    u, v = run(tri, u,v)
    
    with open(tri.filename, 'w', newline='') as file:
        writer = csv.writer(file, delimiter=' ')
        for i in range(n*m):
            writer.writerow([u[i], v[i]])

#
def run_load(tri):
    n = tri.Nx
    m = tri.Ny
    
    u = np.zeros(n*m)
    v = np.zeros(n*m)
        
    with open(tri.filename, newline='') as file:
        reader = csv.reader(file)

        for i in range(n*m):
            ui, vi = next(reader)[0].split(' ')
            u[i] = float(ui)
            v[i] = float(vi)

    u, v = run(tri, u,v)
    
    with open(tri.filename, 'w', newline='') as file:
        writer = csv.writer(file, delimiter=' ')
        for i in range(n*m):
            writer.writerow([u[i], v[i]])
    
    
def run(tri,u,v):
    rhs = make_rhs(tri,u,v)
    M = Dpsi_linOp(tri)

    for i in range(n_trials): 
        print("trial %d of %d"%(i+1, n_trials))
        rhs = make_rhs(tri, u, v)
        stream, exit_flag = bicgstab(M, rhs, atol=stab_tol)
        u, v = uv_approx(tri, u, v, stream)
        
        if i % plot_mod == 0:
            make_plots(tri, u, v)
            
    return u, v
    
# -----------------------------------------------------------------------------
# rhs = c0*A - c1*(B - C)    

# A = u_S - u_N + v_E - v_W
# B = v_C * (u_E + u_W + u_N + u_S)
# C = u_C * (v_E + v_W + v_N + v_S)


def make_rhs(tri, u, v): #
    n = tri.Nx
    m = tri.Ny
    h = tri.h

    yL = tri.yL
    slope = tri.slope
    U = tri.U
    
    c0 = 3*tri.h
    c1 = 0.5 * tri.h**2 * tri.Re
    
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


def uv_approx(tri, u, v, psi):
    
    n = tri.Nx
    m = tri.Ny
    h = tri.h
    
    yL = tri.yL
    slope = tri.slope
    U = tri.U
    
    c2 = 3/(4*h)
    c3 = 1/4
    
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
                psi_E = psi[j*n + i+1] #@ used to be +i
            
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
                v_W = v[j*n + i-1]
                psi_W = psi[j*n + i-1]
   
            new_u[k] = c2 * (psi_N - psi_S) - c3 * (u_N + u_S)
            new_v[k] = -c2 * (psi_E - psi_W) - c3 * (v_E + v_W)
    
    return new_u, new_v

# -----------------------------------------------------------------------------

class Dpsi_linOp(LinearOperator):

    def __init__(self, tri):
        #  globals: Nx, Ny, h, slope, yL
        self.tri = tri
        nm = tri.Nx * tri.Ny
        self.shape = ((nm, nm))        
        self.dtype = np.dtype('f8')
        self.mv = np.zeros(nm) 
        
    def _matvec(self, psi): # v:= psi[i*n + j]  
        n = self.tri.Nx
        m = self.tri.Ny
        h = self.tri.h
        slope = self.tri.slope
        yL = self.tri.yL

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
                
                psi_C = psi[k]
        
                #North (i, j+1)
                y_N = h*(j+1)
                if j+1 == m-1: 
                    psi_N = 0
                else:
                    psi_N = psi[(j+1)*n + i]
                    
                #East (i+1, j)
                x_E = h*(i+1)
                if i+1 == n-1 or y <= slope * x_E - yL:
                    psi_E = 0
                else:
                    psi_E = psi[j*n + (i+1)]
                        
                #South (i, j-1)
                y_S = h*(j-1)
                if j-1 == 0 or y_S <= slope*x - yL or y_S <= -slope*x + yL: 
                    psi_S = 0
                else:
                    psi_S = psi[(j-1)*n + i]
                
                #West (i-1,j)
                x_W = h*(i-1)
                if i-1 == 0 or y <= -slope * x_W - yL: 
                    psi_W = 0
                else:
                    psi_W = psi[j*n + (i-1)]
                    
                #NorthEast (i+1, j+1)
                if i+1 == n-1 or j+1 == m-1 or y_N <= slope * x_E - yL: 
                    psi_NE = 0
                else:
                    psi_NE  = psi[(j+1)*n + (i+1)]
                
                #NorthWest (i-1, j+1)
                if i-1 == 0 or j+1 == m-1 or y_N <= -slope * x_W - yL:
                    psi_NW = 0
                else:
                    psi_NW = psi[(j+1)*n + (i-1)]
                
                #SouthEast (i+1, j-1)
                if i+1 == n-1 or j-1 == 0 or y_S <= slope*x_E - yL or y_S <= -slope*x_E + yL: 
                    psi_SE = 0
                else:
                    psi_SE = psi[(j-1)*n + i+1]
                
                #SouthWest (i-1, j-1)
                if i-1 == 0 or j-1 == 0 or y_S <= slope*x_W - yL or y_S <= -slope*x_W + yL:
                    psi_SW = 0
                else:
                    psi_SW = psi[(j-1)*n + i-1]

            
                self.mv[k] = 28*psi_C - 8*psi_S -8*psi_W -8*psi_E -8*psi_N + psi_SW + psi_SE + psi_NW + psi_NE


        return self.mv

#-----------------------------------------------------------------------------
def load_plot(N):
    tri = biswasKalita(N)
    n = tri.Nx
    m = tri.Ny

    with open(tri.filename, newline='') as file:
        reader = csv.reader(file)
        u = np.ones(n*m)
        v = np.ones(n*m)
        
        for i in range(n*m):
            ui, vi = next(reader)[0].split(' ')
            u[i] = float(ui)
            v[i] = float(vi)
    make_plots(tri, u,v)


def make_plots(tri, u, v):
    n = tri.Nx
    m = tri.Ny

# Grid domain
    xs = np.linspace(tri.x0, tri.xL, n)
    ys = np.linspace(tri.y0, tri.yL, m)

# velocity 
    u_2D = u.reshape((m,n))
    v_2D = v.reshape((m,n))
    
    ax_labels = ['$x$', '$y$']
    title = 'Velocity ($N=%d$)'%(tri.N)
    graph.plot_stream(u_2D, v_2D, xs, ys, title, ax_labels)
    
# vorticity 
    # uy_2D = np.gradient(u_2D, tri.h, axis=0)
    # vx_2D = np.gradient(v_2D, tri.h, axis=1)
    # w = np.zeros((m,n))
    # for j in range(m):
    #     for i in range(n):   
    #         w[j,i] = vx_2D[j,i] - uy_2D[j,i]
    #         if w[j,i] == 0:
    #             w[j,i] = None

    # ax_labels = ['$\omega = v_x - u_y$', '$x$', '$y$']
    # title = 'Vorticity ($N=%d$)'%(tri.N)
    # graph.plot_contour(w, xs, ys, title, ax_labels)

N = 500
mytri = biswasKalita(N)