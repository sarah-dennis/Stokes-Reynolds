#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 22 14:59:46 2024

@author: sarahdennis
"""

import numpy as np

from scipy.interpolate import interpn


import graphics

import stokes_readwrite as rw

# from stokes_solver_tri_spLU import run_spLU
from stokes_solver_BFS_spLU import run_spLU 

write_mod = 250
error_mod = 250

import stokes_examples as examples

# example = examples.biswasEx
# example = examples.zeroReynEx

# example = examples.bfsEx1
# example = examples.bfsEx2
# example = examples.bfsEx3

# example = examples.bfs_biswasLowRe
example = examples.bfs_biswasLowerRe


#------------------------------------------------------------------------------
def new_run(N, iters):
    ex = example(N)

    u_init = np.zeros(ex.Nx * ex.Ny)
    v_init = np.zeros(ex.Nx * ex.Ny)
    psi_init = np.zeros(ex.Nx * ex.Ny)
    past_iters = 0

    u, v, psi = run_spLU(ex, u_init, v_init, psi_init, iters, past_iters, error_mod, write_mod)
    # psi = unmirror_boundary(ex, psi)
    rw.write_solution(ex, u, v, psi, iters)
                                                                                                                                                                                                                                                                             
def load_run(N, iters):
    ex = example(N)

    u, v, psi, past_iters = rw.read_solution(ex.filename+".csv", ex.Nx*ex.Ny)
    
    u, v, psi = run_spLU(ex, u, v, psi, iters, past_iters, error_mod, write_mod)
    
    rw.write_solution(ex, u, v, psi, iters+past_iters)

def load_scale(N_load, N_scale):
    ex_load = example(N_load)
    ex_scale = example(N_scale)
    
    points_load = (ex_load.ys, ex_load.xs)
    
    u_load, v_load, psi_load, past_iters = rw.read_solution(ex_load.filename+".csv", ex_load.Ny*ex_load.Nx)
    u_load_2D = u_load.reshape((ex_load.Ny,ex_load.Nx))
    v_load_2D = v_load.reshape((ex_load.Ny,ex_load.Nx))
    psi_load_2D = psi_load.reshape((ex_load.Ny,ex_load.Nx))


    points_scale = np.meshgrid(ex_scale.ys, ex_scale.xs)
    
    u_scaled_2D = interpn(points_load, u_load_2D, tuple(points_scale), method='cubic')
    v_scaled_2D = interpn(points_load, v_load_2D, tuple(points_scale), method='cubic')
    psi_scaled_2D = interpn(points_load, psi_load_2D, tuple(points_scale), method='cubic')

    u_scaled = u_scaled_2D.ravel(order='F')
    v_scaled = v_scaled_2D.ravel(order='F')
    psi_scaled = psi_scaled_2D.ravel(order='F')

    rw.write_solution(ex_scale, u_scaled, v_scaled, psi_scaled, 0)
    
#------------------------------------------------------------------------------
# PLOTTING 
#------------------------------------------------------------------------------
def load_plot(N):
    ex = example(N)
    u, v, psi, past_iters = rw.read_solution(ex.filename+".csv", ex.Nx * ex.Ny)
    
    n = ex.Nx
    m = ex.Ny
# Grid domain
    xs = ex.xs
    ys = ex.ys

# Stream:
    stream_2D = psi.reshape((ex.Ny,ex.Nx))
    ax_labels = ['$\psi(x,y) = \int u dy + \int v dx$', '$x$', '$y$']
    title = 'Stream ($N=%d$)'%(ex.N)
    log_cmap = True
    graphics.plot_contour_mesh(stream_2D, xs, ys, title, ax_labels, log_cmap)
    
#  Velocity: 
    u_2D = u.reshape((m,n))
    v_2D = v.reshape((m,n))
    ax_labels = ['$|(u,v)|_2$','$x$', '$y$']
    title = 'Velocity ($N=%d$)'%(ex.N)
    ax_labels = ['$\psi(x,y)$ : $u = \psi_y$, $v = -\psi_x$','$x$', '$y$']
    log_cmap = True

    graphics.plot_stream_heat(u_2D, v_2D, xs, ys, stream_2D, title, ax_labels, log_cmap)
    # graphics.plot_quiver_height(u_2D, v_2D, ex.hs, xs, ys, title, ax_labels[1:])
    
#  Vorticity: 
    uy_2D = np.gradient(u_2D, ex.dx, axis=0)
    vx_2D = np.gradient(v_2D, ex.dx, axis=1)
    w = np.zeros((m,n))
    for j in range(m):
        for i in range(n):   
            w[j,i] = vx_2D[j,i] - uy_2D[j,i]
    ax_labels = ['$\omega(x,y) = -( \psi_{xx} + \psi_{yy})$', '$x$', '$y$']
    title = 'Vorticity ($N=%d$)'%(ex.N)
    log_cmap=False
    graphics.plot_contour_mesh(w, xs, ys, title, ax_labels, log_cmap)
    

    
    