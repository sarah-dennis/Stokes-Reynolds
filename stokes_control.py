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

# from stokes_solver_spLU import run_spLU
from stokes_solver_BFS_spLU import run_spLU 

write_mod = 250
error_mod = 250

import stokes_examples as examples
# example = examples.biswasEx
# example = examples.zeroReynEx
# example = examples.bfsEx1
# example = examples.bfsEx2
example = examples.bfsEx3



#------------------------------------------------------------------------------
def new_run(N, iters):
    tri = example(N)

    u_init = np.zeros(tri.Nx * tri.Ny)
    v_init = np.zeros(tri.Nx * tri.Ny)
    psi_init = np.zeros(tri.Nx * tri.Ny)
    past_iters = 0

    u, v, psi = run_spLU(tri, u_init, v_init, psi_init, iters, past_iters, error_mod, write_mod)
    # psi = unmirror_boundary(tri, psi)
    rw.write_solution(tri, u, v, psi, iters)
                                                                                                                                                                                                                                                                             
def load_run(N, iters):
    tri = example(N)

    u, v, psi, past_iters = rw.read_solution(tri.filename+".csv", tri.Nx*tri.Ny)
    
    u, v, psi = run_spLU(tri, u, v, psi, iters, past_iters, error_mod, write_mod)
    
    rw.write_solution(tri, u, v, psi, iters+past_iters)

def load_scale(N_load, N_scale):
    tri_load = example(N_load)
    tri_scale = example(N_scale)
    
    points_load = (tri_load.ys, tri_load.xs)
    
    u_load, v_load, psi_load, past_iters = rw.read_solution(tri_load.filename+".csv", tri_load.Ny*tri_load.Nx)
    u_load_2D = u_load.reshape((tri_load.Ny,tri_load.Nx), order='F')
    v_load_2D = v_load.reshape((tri_load.Ny,tri_load.Nx), order='F')
    psi_load_2D = psi_load.reshape((tri_load.Ny,tri_load.Nx), order='F')


    points_scale = np.meshgrid(tri_scale.ys, tri_scale.xs)
    
    u_scaled_2D = interpn(points_load, u_load_2D, tuple(points_scale), method='linear')
    v_scaled_2D = interpn(points_load, v_load_2D, tuple(points_scale), method='linear')
    psi_scaled_2D = interpn(points_load, psi_load_2D, tuple(points_scale), method='linear')

    u_scaled = u_scaled_2D.ravel()
    v_scaled = v_scaled_2D.ravel()
    psi_scaled = psi_scaled_2D.ravel()

    rw.write_solution(tri_scale, u_scaled, v_scaled, psi_scaled, 0)
    
#------------------------------------------------------------------------------
# PLOTTING 
#------------------------------------------------------------------------------
def load_plot(N):
    tri = example(N)
    u, v, psi, past_iters = rw.read_solution(tri.filename+".csv", tri.Nx * tri.Ny)

    n = tri.Nx
    m = tri.Ny

# Grid domain
    xs = tri.xs
    ys = tri.ys

# Stream: Psi(x,y) heat & contour
    stream_2D = psi.reshape((tri.Ny,tri.Nx))
    ax_labels = ['$\psi(x,y)$ : $u = \psi_y$, $v = -\psi_x$', '$x$', '$y$']
    title = 'Stream ($N=%d$)'%(tri.N)
    graphics.plot_contour_heat(stream_2D, xs, ys, title, ax_labels)
    
    # graphics.plot_2D(stream_2D[:,0], ys, '$\psi(0,y)$', ["$y$", "$\psi$"]) 
    # graphics.plot_2D(stream_2D[:,-1], ys, '$\psi(1,y)$', ["$y$", "$\psi$"]) 
    
#  Velocity: (U, V)  streamplot
    u_2D = u.reshape((m,n))
    v_2D = v.reshape((m,n))
    
    # graphics.plot_2D(u_2D[:,0], ys, '$u(0,y)$', ["$y$", "$u$"])
    # graphics.plot_2D(u_2D[:,-1], ys, '$u(0,y)$', ["$y$", "$u$"])
    
    ax_labels = ['$|(u,v)|_2$','$x$', '$y$']
    title = 'Velocity ($N=%d$)'%(tri.N)
    ax_labels = ['$\psi(x,y)$ : $u = \psi_y$, $v = -\psi_x$','$x$', '$y$']
    graphics.plot_stream_heat(u_2D, v_2D, xs, ys, stream_2D, title, ax_labels)
    graphics.plot_quiver(u_2D, v_2D, xs, ys, stream_2D, title, ax_labels)

#  Vorticity: w = vx - uy heat & contour
    uy_2D = np.gradient(u_2D, tri.dx, axis=0)
    vx_2D = np.gradient(v_2D, tri.dx, axis=1)
    w = np.zeros((m,n))
    for j in range(m):
        for i in range(n):   
            w[j,i] = vx_2D[j,i] - uy_2D[j,i]
    ax_labels = ['$\omega(x,y) = -( \psi_{xx} + \psi_{yy})$', '$x$', '$y$']
    title = 'Vorticity ($N=%d$)'%(tri.N)
    graphics.plot_contour_heat(w, xs, ys, title, ax_labels)
