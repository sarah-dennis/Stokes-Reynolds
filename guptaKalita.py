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

from scipy.sparse.linalg import LinearOperator
from scipy.sparse.linalg import bicgstab

bicgstab_rtol = 1e-10

plot_mod = 50
write_mod = 25
error_mod = 50

class triangle():
    def __init__(self, x0, xL, y0, yL, U, Re, N):
        # N even --> (xL-x0)/2 is triangle vertex
        # slope dividing 2N  --> maximal true boundary points
        self.x0 = x0
        self.xL = xL
        self.y0 = y0
        self.yL = yL
        self.slope = 2*yL/xL
        self.U = U
        self.Re = Re
        self.N = N 
        self.h = 1/N
        self.n = xL*N + 1
        self.m = yL*N + 1
        self.filename = "guptaKalita_N%d.csv"%N

class biswasEx(triangle):
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
    # tri = triangle(x0, xL, y0, yL, U, Re, N)
    tri = biswasEx(N)
    
    n = tri.n
    m = tri.m

    u_init = np.zeros(n*m)
    v_init = np.zeros(n*m)
    psi_init = np.zeros(n*m)
    past_iters = 0

    u, v, psi = run(tri, u_init, v_init, psi_init, iters, past_iters)
    
    write_solution(tri.filename, n*m, u, v, psi, iters)
#
def run_load(N, iters):
    tri = biswasEx(N)
    n = tri.n
    m = tri.m
    
    u, v, psi, past_iters = read_solution(tri.filename, n*m)
    
    u, v, psi = run(tri, u, v, psi, iters, past_iters)
    
    write_solution(tri.filename, n*m, u, v, psi, iters+past_iters)

def plot_load(N):
    tri = biswasEx(N)
    u, v, psi, past_iters = read_solution(tri.filename, tri.n * tri.m)
    make_plots(tri, u, v, psi, past_iters)
    
#------------------------------------------------------------------------------
def run(tri, u, v, past_psi, iters, past_iters):
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 
    M = Dpsi_linOp(tri)
    errs = np.zeros(1 + (iters//error_mod))
    
    for i in range(iters): 
        
        t0 = time.time()
        rhs = update_rhs(tri, u, v)
        
        psi, exit_flag = bicgstab(M, rhs, tol=bicgstab_rtol)
        
        u, v = uv_approx(tri, u, v, psi)
        psi = psi_unmirror_boundary(tri.n, tri.m, tri.slope, psi)
        
        tf = time.time()
        
        if i % error_mod == 0: 
            err_i = np.max(np.abs(psi - past_psi))
            errs[i//error_mod] = err_i
            print("k=%d of %d"%(i+past_iters+1, iters+past_iters))
            print("  time: %.3f s"%(tf-t0))
            print("  error: %.5e psi"%err_i)      
            
        if i % plot_mod == 0:
            make_plots(tri, u, v, psi, i+1 + past_iters)
        
        if i % write_mod == 0:
            write_solution(tri.filename, tri.n*tri.m, u, v, psi, i+1+past_iters)
        
        

        past_psi = psi
    
    trials_max_err = np.max(errs)
    print("max error: %.3e psi"%trials_max_err)
    return u, v, psi

# -----------------------------------------------------------------------------
# rhs = c0*A - c1*(B - C)    

# A = u_S - u_N + v_E - v_W
# B = v_C * (u_E + u_W + u_N + u_S)
# C = u_C * (v_E + v_W + v_N + v_S)

def update_rhs(tri, u, v): #
    
    n = tri.n
    m = tri.m
    
    rhs = np.zeros(n*m)

    h = tri.h
    yL = tri.yL
    slope = tri.slope
    U = tri.U
    
    c0 = 3 * tri.h
    c1 = 0.5 * tri.h**2 * tri.Re
    
    for k in range(n*m):
        i = k % n
        j = k // n
        
        x = h*i
        y = h*j
        
        dy = slope*x
        
        # boundary & dead zones
        if i == 0 or j == 0 or i == n-1 or j == m-1 or y <= -dy + yL or y <= dy - yL:
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
    
    h = tri.h
    yL = tri.yL
    slope = tri.slope
    U = tri.U
    
    c2 = 3/(4*h)
    c3 = 1/4

    for k in range(n*m):
        i = k % n
        j = k // n
        
        x = h * i
        y = h * j
        dy = slope * x
    
        # y=yL boundary
        if j == m-1: 
            u[k] = U
            v[k] = 0 
            
        # other side boundaries & dead zones
        elif i == 0 or j == 0 or i == n-1 or y <= -dy + yL or y <= dy - yL:
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

class Dpsi_linOp(LinearOperator):

    def __init__(self, tri):
        self.tri = tri
        nm = tri.n * tri.m
        self.shape = ((nm, nm))        
        self.dtype = np.dtype('f8')
        self.mv = np.zeros(nm) 
        
    def _matvec(self, psi): # v:= psi[j*n + i]  
        n = self.tri.n
        m = self.tri.m
        
        h = self.tri.h
        slope = self.tri.slope
        yL = self.tri.yL

        # adjust boundary psi's 
        psi = psi_mirror_boundary(n, m, slope, psi)
        
        for k in range(n*m):
            i = k % n
            j = k // n
            
            x = h * i
            y = h * j
            
            dy = slope*x

            # cavity boundary & dead zones
            if i == 0 or i == n-1 or j == 0 or j == m-1 or y <= -dy + yL or y <= dy - yL:
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


def psi_mirror_boundary(n, m, slope, psi):
    i_mid = int((n-1)/2)

    for j in range(m-1):
        
        i_step = int(j//slope) #distance to i_mid \<-|->/

        k_left = j*n + i_mid - i_step #k
        k_right = j*n + i_mid + i_step

        psi[k_left-1] = -psi[k_left+1]
        psi[k_right+1] = -psi[k_right-1]
    
    return psi
        
def psi_unmirror_boundary(n, m, slope, psi):
    i_mid = int((n-1)/2)
    
    for j in range(m-1):
        
        i_step = int(j//slope)
    
        k_left = j*n + i_mid - i_step
        k_right = j*n + i_mid + i_step
    
        psi[k_left-1] = 0
        psi[k_right+1] = 0
    
    return psi

#-----------------------------------------------------------------------------

def make_plots(tri, u, v, stream, iters):
    n = tri.n
    m = tri.m

# Grid domain
    xs = np.linspace(tri.x0, tri.xL, n)
    ys = np.linspace(tri.y0, tri.yL, m)

# Stream: Psi(x,y) heat & contour
    stream_2D = stream.copy().reshape((m,n))
    ax_labels = ['$\psi(x,y)$ : $u = \psi_y$, $v = \psi_x$', '$x$', '$y$']
    title = 'Stream ($N=%d$, $k=%d$)'%(tri.N, iters)
    plot_heat_contour(stream_2D, xs, ys, title, ax_labels, True)

    
# # Velocity: (U, V)  streamplot
#     u_2D = u.copy().reshape((m,n))
#     v_2D = v.copy().reshape((m,n))
    
#     ax_labels = ['$|(u,v)|_2$','$x$', '$y$', ]
#     title = 'Velocity ($N=%d$, $k=%d$)'%(tri.N, iters+1)
#     plot_stream(u_2D, v_2D, xs, ys, title, ax_labels)

# # Vorticity: w = vx - uy heat & contour
#     uy_2D = np.gradient(u_2D, tri.h, axis=0)
#     vx_2D = np.gradient(v_2D, tri.h, axis=1)
#     w = np.zeros((m,n))
#     for j in range(m):
#         for i in range(n):   
#             w[j,i] = vx_2D[j,i] - uy_2D[j,i]
#     ax_labels = ['$\omega(x,y)$', '$x$', '$y$']
#     title = 'Vorticity ($N=%d$, $k=%d$)'%(tri.N, iters)
#     plot_heat_contour(w, xs, ys, title, ax_labels, False)

def plot_heat_contour(zs, xs, ys, title, labels, veriLines):
    pp.rcParams['figure.dpi'] = 500
    pp.figure()
    
    X, Y = np.meshgrid(xs, ys)
    norm_symLog = colors.SymLogNorm(linthresh=bicgstab_rtol, linscale=0.35)
    color_plot = pp.pcolor(X, Y, zs, cmap='Spectral_r', norm=norm_symLog)
    
    pp.colorbar(color_plot, label=labels[0])
    

    n_contours = max(zs.shape)

    pp.rcParams["lines.linewidth"] = .15
    pp.contour(X, Y, zs, n_contours, colors='white')
    
    #line plots
    if veriLines:
        #center line
        pp.plot([0.5, 0.5], [0, 2], '-k')
        
        # vortex centers
        pp.plot([0, 1], [1.802, 1.802], '-k')
        pp.plot([0, 1], [.905, .905], '-k')
        pp.plot([0, 1], [.449, .449], '-k')
        
        # vortex dividers
        pp.plot([0, 1], [1.041, 1.041], '-k')
        pp.plot([0, 1], [.517, .517], '-k')
        pp.plot([0, 1], [.255, .255], '-k')
    
    
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
    magV = np.sqrt(vx**2 + vy**2)
    stream_plot=pp.streamplot(xs, ys, vx, vy, stream_density, linewidth=0.5, color=magV, cmap='Spectral_r', broken_streamlines=False)
    pp.colorbar(stream_plot.lines, label=ax_labels[0])
    
    pp.title(title, fontweight="bold")
    pp.xlabel(ax_labels[1])
    pp.ylabel(ax_labels[2])
    ax = pp.gca()

    ax.set_aspect('equal')
    ax.set_ylim(0)
    pp.show()

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
