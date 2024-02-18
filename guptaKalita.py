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

bicgstab_rtol = 1e-5

plot_mod = 1
write_mod = 5

class cavity():
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

class triangle(cavity):
    def __init__(self, N):
        x0 = 0
        xL = 1
        y0 = 0
        yL = 2

        U = 1
        Re = 1
        super().__init__(x0, xL, y0, yL, U, Re, N)

#------------------------------------------------------------------------------
def run_new(N, iters):
    tri = triangle(N)
    n = tri.Nx
    m = tri.Ny
    
    u_init = np.ones(n*m)
    v_init = np.ones(n*m)
    psi_init = None
    past_iters = 0

    u, v, psi = run(tri, u_init, v_init, iters, past_iters, psi_init)
    
    write_solution(tri.filename, tri.Nx*tri.Ny, u, v, psi, iters)
#
def run_load(N, iters):
    tri = triangle(N)
    n = tri.Nx
    m = tri.Ny
    
    old_u, old_v, old_psi, past_iters = read_solution(tri.filename, n*m)
    
    u, v, psi = run(tri, old_u, old_v, iters, past_iters, old_psi)
    
    write_solution(tri.filename, n*m, u, v, psi, iters+past_iters)

def plot_load(N):
    tri = triangle(N)
    u, v, psi, past_iters = read_solution(tri.filename, tri.Nx * tri.Ny)
    make_plots(tri, u, v, psi, past_iters)
    
#------------------------------------------------------------------------------
def run(tri, u, v, iters, past_iters, past_psi):
    rhs = make_rhs(tri,u,v)                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            
    M = Dpsi_linOp(tri)

    for i in range(iters): 
        print("trial %d of %d"%(+past_iters +i+1, iters+past_iters))
        rhs = make_rhs(tri, u, v)
        psi, exit_flag = bicgstab(M, rhs, tol=bicgstab_rtol)
        u, v = uv_approx(tri, u, v, psi)
        
        if i % plot_mod == 0:
            make_plots(tri, u, v, psi, i+past_iters)
        
        if i % write_mod == write_mod - 1:
            write_solution(tri.filename, tri.Nx*tri.Ny, u, v, psi, i+1+past_iters)
        
        
        if past_psi.any(): 
            infNormErr = np.max(np.abs(psi - past_psi))
            print(infNormErr)
            
        past_psi = psi
    return u, v, psi

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
            if i-1 == 0 or y <= -slope * x_W + yL: 
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
                psi_E = psi[j*n + i+1] 
            
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
            if i-1 == 0 or y <= -slope * x_W + yL:
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
                if i-1 == 0 or y <= -slope * x_W + yL:
                    psi_W = 0
                else:
                    psi_W = psi[j*n + (i-1)]
                    
                #NorthEast (i+1, j+1)
                if i+1 == n-1 or j+1 == m-1 or y_N <= slope * x_E - yL: 
                    psi_NE = 0
                else:
                    psi_NE  = psi[(j+1)*n + (i+1)]
                
                #NorthWest (i-1, j+1)
                if i-1 == 0 or j+1 == m-1 or y_N <= -slope * x_W + yL: 
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

def make_plots(tri, u, v, stream, iters):
    n = tri.Nx
    m = tri.Ny

# Grid domain
    xs = np.linspace(tri.x0, tri.xL, n)
    ys = np.linspace(tri.y0, tri.yL, m)

# velocity 
    # u_2D = u.reshape((m,n))
    # v_2D = v.reshape((m,n))
    
    # ax_labels = ['$x$', '$y$']
    # title = 'Velocity ($N=%d$, $k=%d$)'%(tri.N, iters+1)
    # graph.plot_stream(u_2D, v_2D, xs, ys, title, ax_labels)
    
# stream 
    stream_2D = u.reshape((m,n))
    for j in range(m):
        for i in range(n):   
            if stream_2D[j,i] == 0:
                stream_2D[j,i] = None
    ax_labels = ['$\psi(x,y)$', '$x$', '$y$']
    title = 'Stream ($N=%d$, $k=%d$)'%(tri.N, iters+1)
    graph.plot_heat_contour(stream_2D, xs, ys, title, ax_labels)
    # graph.plot_contour(stream_2D, xs, ys, title, ax_labels)
    
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
    # title = 'Vorticity ($N=%d$, $k=%d$)'%(tri.N, iters+1)
    # graph.plot_heat_contour(w, xs, ys, title, ax_labels)

#------------------------------------------------------------------------------
def read_solution(filename, nm):
    u = np.zeros(nm)
    v = np.zeros(nm)
    psi = np.zeros(nm)
    with open(filename, newline='') as file:
        reader = csv.reader(file)

        for i in range(nm):
            ui, vi, psii = next(reader)[0].split(' ')
            u[i] = float(ui)
            v[i] = float(vi)
            psi[i] = float(psii)
        past_iters = int(next(reader)[0])
        file.close()
    return u, v, psi, past_iters

def write_solution(filename, nm, u, v, psi, iters):

    with open(filename, 'w', newline='') as file:
        writer = csv.writer(file, delimiter=' ')
        for i in range(nm):
            writer.writerow([u[i], v[i], psi[i]])
        
        writer.writerow([iters])
        print("Saved to csv: k=%d"%iters)
        file.close()
