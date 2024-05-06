#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 22 14:59:46 2024

@author: sarahdennis
"""
import csv
import time
import numpy as np

from matplotlib import pyplot as pp
from matplotlib import colors
from matplotlib import patches

from scipy.sparse.linalg import LinearOperator
from scipy.sparse import csc_matrix
from scipy.sparse.linalg import bicgstab
from scipy.linalg import lu_factor, lu_solve
from scipy.sparse.linalg import splu

from scipy.signal import argrelextrema as relEx

# from scipy.ndimage import zoom

from scipy.interpolate import interpn
# from scipy.interpolate import griddata




bicgstab_rtol = 1e-7

bicgstab_rtol = 1e-6

plot_mod = 20
write_mod = 20
error_mod = 20


class triangle():
    def __init__(self, x0, xL, y0, yL, U, Re, N, filestr):
        # N even --> (xL-x0)/2 is triangle vertex
        # slope dividing 2N  --> maximal true boundary points
        self.x0 = x0
        self.xL = xL
        self.y0 = y0
        self.yL = yL
        self.slope = int(2*yL/xL)
        self.U = U
        self.Re = Re
        self.N = N 
        self.h = 1/N
        self.n = xL*N + 1
        self.m = yL*N + 1
        self.filename = filestr
        self.xs = np.linspace(x0, xL, self.n)
        self.ys = np.linspace(y0, yL, self.m)

class biswasEx(triangle):
    def __init__(self, N):
        x0 = 0
        xL = 1
        y0 = 0
        yL = 2
        U = 1
        Re = 1
        filestr = "stokes_N%d.csv"%(N)
        super().__init__(x0, xL, y0, yL, U, Re, N, filestr)
        
#------------------------------------------------------------------------------
def run_new(N, iters):
    # tri = triangle(x0, xL, y0, yL, U, Re, N)
    tri = biswasEx(N)

    nm = tri.n * tri.m

    u_init = np.zeros(nm)
    v_init = np.zeros(nm)
    psi_init = np.zeros(nm)
    past_iters = 0

    # u, v, psi = run(tri, u_init, v_init, psi_init, iters, past_iters)
    # u, v, psi = run_LU(tri, u_init, v_init, psi_init, iters, past_iters)
    u, v, psi = run_spLU(tri, u_init, v_init, psi_init, iters, past_iters)
    
    psi = psi_unmirror_boundary(tri, psi)
    write_solution(tri.filename+"_LU", nm, u, v, psi, iters)
                                                                                                                                                                                                                                                                             
def run_load(N, iters):
    # tri = triangle(x0, xL, y0, yL, U, Re, N)
    tri = biswasEx(N)

    u, v, psi, past_iters = read_solution(tri.filename, tri.n*tri.m)
    
    # u, v, psi = run(tri, u, v, psi, iters, past_iters)
    # u, v, psi = run_LU(tri, u, v, psi, iters, past_iters)
    u, v, psi = run_spLU(tri, u, v, psi, iters, past_iters)
    
    psi = psi_unmirror_boundary(tri, psi)
    write_solution(tri.filename+"_LU", tri.n*tri.m, u, v, psi, iters+past_iters)


def load_scale(N_load, N_new):
    # tri_load = triangle(x0, xL, y0, yL, U, Re, N)
    tri_load = biswasEx(N_load) #slope 4 example
    
    points_load = (tri_load.ys, tri_load.xs)
    
    u_load, v_load, psi_load, past_iters = read_solution(tri_load.filename, tri_load.m*tri_load.n)
    u_load_2D = u_load.reshape((tri_load.m,tri_load.n), order='F')
    v_load_2D = v_load.reshape((tri_load.m,tri_load.n), order='F')
    psi_load_2D = psi_load.reshape((tri_load.m,tri_load.n), order='F')


    tri_scale = biswasEx(N_new)
    
    # previously...
    # new_shape = (tri_scale.m/tri_load.m, tri_scale.n/tri_load.n)
    # u_scaled_2D = zoom(u_load_2D, new_shape)
    # v_scaled_2D = zoom(v_load_2D, new_shape)
    # psi_scaled_2D = zoom(psi_load_2D, new_shape)
    
    # and now...

    points_scale = np.meshgrid(tri_scale.ys, tri_scale.xs)
    u_scaled_2D = interpn(points_load, u_load_2D, tuple(points_scale), method='linear')
    v_scaled_2D = interpn(points_load, v_load_2D, tuple(points_scale), method='linear')
    psi_scaled_2D = interpn(points_load, psi_load_2D, tuple(points_scale), method='linear')



    u_scaled = u_scaled_2D.ravel()
    v_scaled = v_scaled_2D.ravel()
    psi_scaled = psi_scaled_2D.ravel()

    write_solution(tri_scale.filename, tri_scale.m*tri_scale.n, u_scaled, v_scaled, psi_scaled, 0)
    plot_load(N_new)
    
def plot_load(N):
    # tri = triangle(x0, xL, y0, yL, U, Re, N)
    tri = biswasEx(N)
    u, v, psi, past_iters = read_solution(tri.filename, tri.n * tri.m)
    make_plots(tri, u, v, psi, past_iters)

#------------------------------------------------------------------------------
def read_solution(filename, nm):
    u = np.zeros(nm)
    v = np.zeros(nm)
    psi = np.zeros(nm)
    with open(filename, newline='') as file:
        reader = csv.reader(file)

        for i in range(nm):
            line = next(reader)
            ui, vi, psii = line[0].split(' ')
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
        print("  saved csv")
        file.close()    

#------------------------------------------------------------------------------
def run(tri, u, v, past_psi, iters, past_iters):
    ##TODO: speed-up ?
    
    # precondition M since it never changes
    # bicgstab accepts a pseudo inverse for M as preconditioning
    
    M = Dpsi_linOp(tri)
    
    for i in range(iters): 
        
        t0 = time.time()
        rhs = update_rhs(tri, u, v)
        
        psi, exit_flag = bicgstab(M, rhs, tol=bicgstab_rtol)
        
        
        u, v = uv_approx(tri, u, v, psi)
        
        psi = psi_mirror_boundary(tri, psi)
        
        tf = time.time()
        if i % error_mod == 0: 
            past_psi = psi_unmirror_boundary(tri, past_psi)
            psi = psi_unmirror_boundary(tri, psi)
            
            err_i = np.max(np.abs(psi - past_psi))
            
            print("k=%d of %d"%(i+past_iters+1, iters+past_iters))
            print("  time: %.3f s"%(tf-t0))
            print("  error: %.5e psi"%err_i)
            
            if err_i < bicgstab_rtol:
                print('warning: lower bicgstab_rtol')
            
        if i % plot_mod == 0:
            psi = psi_unmirror_boundary(tri, psi)
            make_plots(tri, u, v, psi, i+1 + past_iters)
        
        if i % write_mod == 0:
            psi = psi_unmirror_boundary(tri, psi)
            write_solution(tri.filename, tri.n*tri.m, u, v, psi, i+1+past_iters)


        past_psi = psi

    return u, v, psi

def run_LU(tri, u, v, past_psi, iters, past_iters):
    
    M = Dpsi_matrixBuild(tri)
    LU, piv = lu_factor(M)
    for i in range(iters): 
        
        t0 = time.time()
        rhs = update_rhs(tri, u, v)

        psi = lu_solve((LU, piv), rhs)
        
        u, v = uv_approx(tri, u, v, psi)
        
        psi = psi_mirror_boundary(tri, psi)
        
        tf = time.time()
        
        if i % error_mod == 0: 
            past_psi = psi_unmirror_boundary(tri, past_psi)
            psi = psi_unmirror_boundary(tri, psi)
            
            err_i = np.max(np.abs(psi - past_psi))
            
            print("k=%d of %d"%(i+past_iters+1, iters+past_iters))
            print("  time: %.3f s"%(tf-t0))
            print("  error: %.5e psi"%err_i)
            
        if i % plot_mod == 0:
            psi = psi_unmirror_boundary(tri, psi)
            make_plots(tri, u, v, psi, i+1 + past_iters)
        
        if i % write_mod == 0:
            psi = psi_unmirror_boundary(tri, psi)
            write_solution(tri.filename+"_LU", tri.n*tri.m, u, v, psi, i+1+past_iters)


        past_psi = psi

    return u, v, psi

def run_spLU(tri, u, v, past_psi, iters, past_iters):
    
    M = Dpsi_cscmatrixBuild(tri)
    LU = splu(M)
    
    for i in range(iters): 
        
        t0 = time.time()
        rhs = update_rhs(tri, u, v)

        psi = LU.solve(rhs)
        
        u, v = uv_approx(tri, u, v, psi)
        
        psi = psi_mirror_boundary(tri, psi)
        
        tf = time.time()
        
        if i % error_mod == 0: 
            past_psi = psi_unmirror_boundary(tri, past_psi)
            psi = psi_unmirror_boundary(tri, psi)
            
            err_i = np.max(np.abs(psi - past_psi))
            
            print("k=%d of %d"%(i+past_iters+1, iters+past_iters))
            print("  time: %.3f s"%(tf-t0))
            print("  error: %.5e psi"%err_i)
            
        if i % plot_mod == 0:
            psi = psi_unmirror_boundary(tri, psi)
            make_plots(tri, u, v, psi, i+1 + past_iters)
        
        if i % write_mod == 0:
            psi = psi_unmirror_boundary(tri, psi)
            write_solution(tri.filename+"_LU", tri.n*tri.m, u, v, psi, i+1+past_iters)


        past_psi = psi

    return u, v, psi



# -----------------------------------------------------------------------------
# rhs = c0*A - c1*(B - C)

# A = u_S - u_N + v_E - v_W
# B = v_C * (u_E + u_W + u_N + u_S)
# C = u_C * (v_E + v_W + v_N + v_S)

def update_rhs(tri, u, v): #
    
    n = tri.n 
    nc = n//2 # i-index of triangle apex
    m = tri.m
    
    rhs = np.zeros(n*m)

    # h = tri.h
    # yL = tri.yL
    slope = tri.slope
    U = tri.U
    
    c0 = 3 * tri.h
    c1 = 0.5 * tri.h**2 * tri.Re
    
    for k in range(n*m):
        i = k % n
        j = k // n

        # boundary & dead zones
        if i == 0 or j == 0 or i == n-1 or j == m-1:
            rhs[k] = 0
               
        elif (i < nc and j < (m-1) - slope*i) or (i > nc and j < slope * (i - nc)):
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
                k_N = (j+1)*n + i
                u_N = u[k_N]
                v_N = v[k_N]
                
            # East (i+1, j)
            if i+1 == n-1: 
                u_E = 0
                v_E = 0
            else:
                k_E = j*n + i + 1
                u_E = u[k_E]
                v_E = v[k_E]
            
            # South (i, j-1)
            if j-1 == 0: 
                u_S = 0
                v_S = 0
            else:
                k_S = (j-1)*n + i
                u_S = u[k_S]
                v_S = v[k_S]

                
            # West (i-1, j)  
            if i-1 == 0: 
                u_W = 0
                v_W = 0
            else:
                k_W = j*n + i - 1
                u_W = u[k_W]
                v_W = v[k_W]
             
            A = u_S - u_N + v_E - v_W
            B = v_C * (u_E + u_W + u_N + u_S)
            C = u_C * (v_E + v_W + v_N + v_S)

            rhs[k] = c0 * A + c1 * (B - C)
            
    return rhs

# -----------------------------------------------------------------------------

# u[k] = c2 * (psi_N - psi_S) - c3 * (u_N + u_S)
# v[k] = -c2 * (psi_E - psi_W) - c3 * (v_E + v_W)

def uv_approx(tri, u, v, psi):
    
    n = tri.n
    m = tri.m

    nc = n//2
    h = tri.h
    slope = tri.slope
    U = tri.U
    
    c2 = 3/(4*h)
    c3 = 1/4

    for k in range(n*m):
        i = k % n
        j = k // n

        # y=yL moving boundary
        if j == m-1: 
            u[k] = U
            v[k] = 0 
            
        # other boundaries & dead zones
        elif j == 0 or i == 0 or i == n-1:
            u[k] = 0
            v[k] = 0 
        
        elif (i < nc and j < (m-1) - slope*i) or (i > nc and j < slope * (i - nc)):
            u[k] = 0
            v[k] = 0 
                
        else: #interior
            # (u,v) at 4 point stencil

            # North (i, j+1)
            if j+1 == m-1: 
                u_N = U
                psi_N = 0
            else:
                k_N = (j+1)*n + i
                u_N = u[k_N]
                psi_N = psi[k_N]
                
            # East (i+1, j)
            if i+1 == n-1 : 
                v_E = 0
                psi_E = 0 
            else:
                k_E = j*n + i + 1
                v_E = v[k_E]
                psi_E = psi[k_E] 
            
            # South (i, j-1)
            if j-1 == 0 : 
                u_S = 0
                psi_S = 0 
            else:
                k_S = (j-1)*n + i
                u_S = u[k_S]
                psi_S = psi[k_S]
                
            # West (i-1, j)  
            if i-1 == 0 :
                v_W = 0
                psi_W = 0 
            else:
                k_W = j*n + i - 1
                v_W = v[k_W]
                psi_W = psi[k_W]
   
            u[k] = c2 * (psi_N - psi_S) - c3 * (u_N + u_S)
            v[k] = -c2 * (psi_E - psi_W) - c3 * (v_E + v_W)
    
    return u, v

# -----------------------------------------------------------------------------
# LinOp dPsi : 9 point discr. 
# -----------------------------------------------------------------------------

class Dpsi_linOp(LinearOperator):

    def __init__(self, tri):
        self.tri = tri
        nm = tri.n * tri.m
        self.shape = ((nm, nm))        
        self.dtype = np.dtype('f8')
        self.mv = np.zeros(nm) 
        
    def _matvec(self, psi): # v:= psi[j*n + i]  
        n = self.tri.n
        nc = n//2
        m = self.tri.m
        
        slope = self.tri.slope

        for k in range(n*m):
            i = k % n
            j = k // n
            
            # cavity boundary & dead zones
            if i == 0 or i == n-1 or j == 0 or j == m-1: 
                self.mv[k] = 0
            
            elif  (i < nc and j < (m-1) - slope*i) or (i > nc and j < slope * (i - nc)):
                self.mv[k] = 0
       
            #interior   
            else: 
                # psi[k] at 9 point stencil         
                
                psi_C = psi[k]
        
                #North (i, j+1)
                if j+1 == m-1: 
                    psi_N = 0
                else:
                    psi_N = psi[(j+1)*n + i]
                    
    
                #South (i, j-1)
                if j-1 == 0 : 
                    psi_S = 0
                else:
                    psi_S = psi[(j-1)*n + i]
                    
                #East (i+1, j)
                if i+1 == n-1:
                    psi_E = 0
                else:
                    psi_E = psi[j*n + i+1]
                    
                #West (i-1,j)
                if i-1 == 0 :
                    psi_W = 0
                else:
                    psi_W = psi[j*n + i-1]
                    
                #NorthEast (i+1, j+1)
                if i+1 == n-1 or j+1 == m-1 : 
                    psi_NE = 0
                else:
                    psi_NE  = psi[(j+1)*n + i+1]
                
                #NorthWest (i-1, j+1)
                if i-1 == 0 or j+1 == m-1: 
                    psi_NW = 0
                else:
                    psi_NW = psi[(j+1)*n + i-1]
                
                #SouthEast (i+1, j-1)
                if i+1 == n-1 or j-1 == 0: 
                    psi_SE = 0
                else:
                    psi_SE = psi[(j-1)*n + i+1]
                
                #SouthWest (i-1, j-1)
                if i-1 == 0 or j-1 == 0:
                    psi_SW = 0
                else:
                    psi_SW = psi[(j-1)*n + i-1]

            
                self.mv[k] = 28*psi_C - 8*(psi_N + psi_S + psi_E + psi_W) + psi_NE + psi_SE + psi_NW + psi_SW

        return self.mv

def Dpsi_matrixBuild(tri):
    m = tri.m
    n = tri.n
    mat = np.zeros((m*n, m*n))
    slope = tri.slope
    nc = n//2
    #fill mat row by row
    for k in range(m*n):
        
    
        i = k % n
        j = k // n
        
        # k = j*n + i
        
        if i == 0 or j == 0 or i == n-1 or j == m-1:
            mat[k][k] = 1
        
        elif  (i < nc and j < (m-1) - slope*i) or (i > nc and j < slope * (i - nc)):
            mat[k][k] = 1

        else: 
            
            mat[k][k] = 28
            
            mat[k][j*n + i - 1] = -8
            mat[k][j*n + i + 1] = -8
            mat[k][(j-1)*n + i] = -8
            mat[k][(j+1)*n + i] = -8
            
            mat[k][(j-1)*n + i-1] = 1
            mat[k][(j-1)*n + i+1] = 1
            mat[k][(j+1)*n + i-1] = 1
            mat[k][(j+1)*n + i+1] = 1
    
    return mat

def Dpsi_cscmatrixBuild(tri):
    m = tri.m
    n = tri.n
    mat = csc_matrix((m*n, m*n))
    slope = tri.slope
    nc = n//2
    #fill mat row by row
    for k in range(m*n):
        
    
        i = k % n
        j = k // n
        
        # k = j*n + i
        
        if i == 0 or j == 0 or i == n-1 or j == m-1:
            mat[k][k] = 1
        
        elif  (i < nc and j < (m-1) - slope*i) or (i > nc and j < slope * (i - nc)):
            mat[k][k] = 1

        else: 
            
            mat[k][k] = 28
            
            mat[k][j*n + i - 1] = -8
            mat[k][j*n + i + 1] = -8
            mat[k][(j-1)*n + i] = -8
            mat[k][(j+1)*n + i] = -8
            
            mat[k][(j-1)*n + i-1] = 1
            mat[k][(j-1)*n + i+1] = 1
            mat[k][(j+1)*n + i-1] = 1
            mat[k][(j+1)*n + i+1] = 1
    
    return mat


# -----------------------------------------------------------------------------
# boundary mirroring

def psi_mirror_boundary(tri, psi):
    n = tri.n 
    m = tri.m 
    slope = tri.slope
    dx = tri.h
    
    i_mid = n//2
    
    for j in range(m-1):
        
        dj = j % slope
        di = int(j//slope)
        
        k_left = j*n + i_mid - di
        k_right = j*n + i_mid + di
        
        if dj == 0: #true boundary point
            psi[k_left] = 0 
            psi[k_right] = 0
        
        else:
    
            psi[k_left-1] = (1 - dj*dx/slope)*psi[k_left]
            psi[k_right+1] = (1 + dj*dx/slope)*psi[k_right] 
             
    return psi
        
def psi_unmirror_boundary(tri, psi):
    n = tri.n
    m = tri.m
    slope = tri.slope
    
    i_mid = n//2

    for j in range(m-1):
        
        di = int(j//slope) 
    
        k_left = j*n + i_mid - di
        k_right = j*n + i_mid + di
        
        psi[k_left-1] = 0
        psi[k_right+1] = 0
    
    return psi



#------------------------------------------------------------------------------
# PLOTTING 
#------------------------------------------------------------------------------
def make_plots(tri, u, v, stream, iters):
    n = tri.n
    m = tri.m

# Grid domain
    xs = tri.xs
    ys = tri.ys

# Stream: Psi(x,y) heat & contour
    stream_2D = stream.copy().reshape((m,n))
    ax_labels = ['$\psi(x,y)$ : $u = \psi_y$, $v = \psi_x$', '$x$', '$y$']
    title = 'Stream ($N=%d$, $k=%d$)'%(tri.N, iters)
    plot_contour_heat(stream_2D, xs, ys, title, ax_labels)
      
#  Velocity: (U, V)  streamplot
    u_2D = u.copy().reshape((m,n))
    v_2D = v.copy().reshape((m,n))
    
    ax_labels = ['$|(u,v)|_2$','$x$', '$y$']
    title = 'Velocity ($N=%d$, $k=%d$)'%(tri.N, iters+1)
    # plot_stream(u_2D, v_2D, xs, ys, title, ax_labels)
    ax_labels = ['$\psi(x,y)$ : $u = \psi_y$, $v = \psi_x$','$x$', '$y$']
    plot_stream_heat(u_2D, v_2D,  xs, ys, stream_2D, title, ax_labels)

#  Vorticity: w = vx - uy heat & contour
    uy_2D = np.gradient(u_2D, tri.h, axis=0)
    vx_2D = np.gradient(v_2D, tri.h, axis=1)
    w = np.zeros((m,n))
    for j in range(m):
        for i in range(n):   
            w[j,i] = vx_2D[j,i] - uy_2D[j,i]
    ax_labels = ['$\omega(x,y) = -( \psi_{xx} + \psi_{yy})$', '$x$', '$y$']
    title = 'Vorticity ($N=%d$, $k=%d$)'%(tri.N, iters)
    plot_contour_heat(w, xs, ys, title, ax_labels)

def plot_contour_heat(zs, xs, ys, title, labels):
    pp.rcParams['figure.dpi'] = 500
    pp.figure()
    
    X, Y = np.meshgrid(xs, ys)

    norm_symLog = colors.SymLogNorm(linthresh=1e-12, linscale=0.35)

    color_plot = pp.pcolor(X, Y, zs, cmap='Spectral_r', norm=norm_symLog)
    
    pp.colorbar(color_plot, label=labels[0])
    

    n_contours = 10

    pp.rcParams["lines.linewidth"] = .15
    pp.contour(X, Y, zs, n_contours, colors='white')
    
    pp.title(title, fontweight="bold")
    pp.xlabel(labels[1])
    pp.ylabel(labels[2])
    
    ax = pp.gca()
    ax.set_aspect('equal', 'box')
    pp.show()

def plot_stream(vx, vy, xs, ys, title, ax_labels):
    
    pp.rcParams['figure.dpi'] = 500
    pp.figure()
    
    X, Y = np.meshgrid(xs, ys)
    
    stream_density=[1,2] #len(ys) = 2 len(xs)
    
    pp.streamplot(xs, ys, vx, vy, stream_density, color='k', linewidth=0.5, broken_streamlines=False, )
    
    pp.title(title, fontweight="bold")
    pp.xlabel(ax_labels[1])
    pp.ylabel(ax_labels[2])
    ax = pp.gca()

    ax.set_aspect('equal')
    ax.set_ylim(0)
    pp.show()
    
def plot_stream_heat(vx, vy, xs, ys, psi, title, ax_labels):
    
    pp.rcParams['figure.dpi'] = 500
    pp.figure()


    X, Y = np.meshgrid(xs, ys)
    
    stream_density=[1,2] #len(ys) = 2 len(xs)
    
    norm_symLog = colors.SymLogNorm(linthresh=1e-12, linscale=0.35)

    stream_plot=pp.streamplot(xs, ys, vx, vy, stream_density, linewidth=0.5, color=psi, cmap='Spectral_r',  norm=norm_symLog, broken_streamlines=False)
    
    ax = pp.gca()
    for art in ax.get_children():
        if not isinstance(art, patches.FancyArrowPatch):
            continue
        art.remove()        
    
    
    
    pp.colorbar(stream_plot.lines, label=ax_labels[0])
    # stream_plot.arrows.FancyArrowPatch.set_arrowstyle("fancy", visible=False)
    
    # magV = np.sqrt(vx**2 + vy**2)
    # stream_plot=pp.streamplot(xs, ys, vx, vy, stream_density, linewidth=0.5, color=magV, cmap='Spectral_r',  norm=norm_symLog, broken_streamlines=False)
    # pp.colorbar(stream_plot.lines, label=ax_labels[0])
    
    pp.title(title, fontweight="bold")
    pp.xlabel(ax_labels[1])
    pp.ylabel(ax_labels[2])

    ax.set_aspect('equal')
    ax.set_ylim(0,2) #min max ys
    pp.show()



#------------------------------------------------------------------------------
# Critiacal points 
#------------------------------------------------------------------------------

def get_boundary(tri, psi):
    n = tri.n
    m = tri.m
    
    i_mid = n//2
    
    left = np.zeros((m,3)) # x, y, psi
    right = np.zeros((m,3)) 
    
    
    for j in range(m):
        
        y = tri.y0 + j*tri.h
        
        dj = j % tri.slope
        di = int(j//tri.slope) 
        
        if dj == 0: #true boundary points
            k_left = j*n + i_mid - di 
            k_right = j*n + i_mid + di
        else: #interior boundary points
            k_left = j*n + i_mid - di + 1  
            k_right = j*n + i_mid + di - 1
        
        x_left = (i_mid - di + 1)*tri.h + tri.x0
        x_right = tri.xL - (i_mid - di - 1)*tri.h
        
        left[j] = [x_left, y, psi[k_left]]
        right[j] = [x_right, y, psi[k_right]]
        
    return left, right

def get_center(tri, psi):
    n = tri.n 
    m = tri.m 
    
    i_mid = n//2
    x = i_mid * tri.h + tri.x0
    
    center = np.zeros((m,3)) # x, y, psi
    
    for j in range(m):
        
        y = tri.y0 + j*tri.h
        
        k = j*n + i_mid
        
        center[j] = [x, y, psi[k]]
        
    return center

def write_criticals(N):
    tri = biswasEx(N)
    u, v, psi, past_iters = read_solution(tri.filename, tri.n * tri.m)
    
    left, right = get_boundary(tri, psi)
    
    crits_filename = 'crits_' + tri.filename 
    with open(crits_filename, 'w', newline='') as file:
        writer = csv.writer(file, delimiter=' ')
    
        # print("Psi sign changes...")
        writer.writerow('xy')
        sign_ref = 0
        for (x, y, p) in left:
            sign_new = np.sign(p)
            if sign_new != 0 and sign_new != sign_ref:
                sign_ref = sign_new
                # print("(x:%.5f, y:%.6f)"%(x,y))
                writer.writerow([x,y])
                
        writer.writerow('xy')
        sign_ref = 0
        for (x, y, p) in right:
            sign_new = np.sign(p)
            if sign_new != 0 and sign_new != sign_ref:
                sign_ref = sign_new
                # print("(x:%.5f, y:%.6f)"%(x,y))
                writer.writerow([x,y])

    
        # print("Psi center-line extrema...")
        writer.writerow('xyp')
        center = get_center(tri, psi)
        max_inds = relEx(center[:,2], np.greater)[0]
        min_inds = relEx(center[:,2], np.less)[0]
        
    
        for i in max_inds:
            x,y,p = center[i]
            # print("(x:%.1f, y:%.6f) p=%.5e"% (x,y,p))
            writer.writerow([x,y,p])
        
        writer.writerow('xyp')
        for i in min_inds:
            x,y,p = center[i]
            # print("(x:%.1f, y:%.6f) p=%.5e"% (x,y,p))
            writer.writerow([x,y,p])

    make_plots(tri, u, v, psi, past_iters)

# def read_criticals(N, k_hs, k_extr):
#     tri = biswasEx(N)
#     crits_filename = 'crits_' + tri.filename 
    
    
#     with open(crits_filename, newline='') as file:
#         reader = csv.reader(file)
#         for line in reader:
#             if line[0] !== 'x':
                
            

# def write_criticals_error(N_a, N_b): # compare N_a grid with N_b grid


    