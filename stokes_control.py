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
import stokes_convergence as cnvg

import stokes_examples as examples

from stokes_solver_tri import run_spLU
# example = examples.tri_Re1
example = examples.tri_Re0

# from stokes_solver_BFS import run_spLU 
# example = examples.bfs_Re10neg4
# example = examples.bfs_Re0

write_mod = 500
error_mod = 500

#------------------------------------------------------------------------------
def new_run(N, iters):
    
    ex = example(N)
    u_init = np.ones(ex.Nx * ex.Ny)
    v_init = np.ones(ex.Nx * ex.Ny)
    psi_init = np.ones(ex.Nx * ex.Ny)
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
def vorticity(ex, u_2D, v_2D):
    uy_2D = np.gradient(u_2D, ex.dx, axis=0)
    vx_2D = np.gradient(v_2D, ex.dx, axis=1)
    w_2D = np.zeros((ex.Ny,ex.Nx))
    for j in range(ex.Ny):
        for i in range(ex.Nx):   
            w_2D[j,i] = vx_2D[j,i] - uy_2D[j,i]
    return w_2D
#------------------------------------------------------------------------------
# PLOTTING 
#------------------------------------------------------------------------------
def load_plot(N):
    ex = example(N)
    u, v, psi, past_iters = rw.read_solution(ex.filename+".csv", ex.Nx * ex.Ny)

# Grid domain
    xs = ex.xs
    ys = ex.ys
    
    stream_2D = psi.reshape((ex.Ny,ex.Nx))
    u_2D = u.reshape((ex.Ny,ex.Nx))
    v_2D = v.reshape((ex.Ny,ex.Nx))

# zoom domain
    x_start_A = 0.3
    x_stop_A = 0.7
    y_start_A = 0
    y_stop_A = .6

    xs_zoom, ys_zoom = grid_zoom_1D(xs, ys, ex, x_start_A, x_stop_A, y_start_A, y_stop_A)
    stream_2D_zoom = grid_zoom_2D(stream_2D, ex, x_start_A, x_stop_A, y_start_A, y_stop_A)
    u_2D_zoom = grid_zoom_2D(u_2D, ex, x_start_A, x_stop_A, y_start_A, y_stop_A)
    v_2D_zoom = grid_zoom_2D(v_2D, ex, x_start_A, x_stop_A, y_start_A, y_stop_A)

# Stream plot:

    ax_labels = ['$\psi(x,y) = \int u dy + \int v dx$', '$x$', '$y$']
    title = 'Stream ($N=%d$, Re$=%.5f$)'%(ex.N, ex.Re)
    graphics.plot_contour_mesh(stream_2D, xs, ys, title, ax_labels, True)
    graphics.plot_contour_mesh(stream_2D_zoom, xs_zoom, ys_zoom, title, ax_labels, True)
    
#  Velocity plot: 

    ax_labels = ['$|(u,v)|_2$','$x$', '$y$']
    title = 'Velocity ($N=%d$, Re$=%.5f$)'%(ex.N, ex.Re)
    ax_labels = ['$\psi(x,y)$ : $u = \psi_y$, $v = -\psi_x$','$x$', '$y$']
    graphics.plot_stream_heat(u_2D, v_2D, xs, ys, stream_2D, title, ax_labels, True) 

    graphics.plot_stream_heat(u_2D_zoom, v_2D_zoom, xs_zoom, ys_zoom, stream_2D_zoom, title, ax_labels, True)
    
#  Vorticity plot: 
    w = vorticity(ex, u_2D, v_2D)

    ax_labels = ['$\omega(x,y) = -( \psi_{xx} + \psi_{yy})$', '$x$', '$y$']
    title = 'Vorticity ($N=%d$, Re$=%.5f$)'%(ex.N, ex.Re)
    graphics.plot_contour(w, xs, ys, title, ax_labels)
                
    x_start_B = 0
    x_stop_B = 0.1
    y_start_B = 1.75
    y_stop_B = 2
    w_zoom = grid_zoom_2D(w, ex, x_start_B, x_stop_B, y_start_B, y_stop_B)     
    xs_zoom, ys_zoom = grid_zoom_1D(xs, ys, ex, x_start_B, x_stop_B, y_start_B, y_stop_B)
    
    graphics.plot_contour(w_zoom, xs_zoom, ys_zoom, title, ax_labels)

def grid_zoom_2D(grid, ex, x_start, x_stop, y_start, y_stop):
    i_0 = int((x_start - ex.x0)/ex.dx)
    i_f = int((x_stop - ex.x0)/ex.dx)
    j_0 = int((y_start - ex.y0)/ex.dy)
    j_f = int((y_stop - ex.y0)/ex.dy)
    return grid[j_0:j_f,i_0:i_f]

def grid_zoom_1D(grid_x, grid_y, ex, x_start, x_stop, y_start, y_stop):
    i_0 = int((x_start - ex.x0)/ex.dx)
    i_f = int((x_stop - ex.x0)/ex.dx)
    j_0 = int((y_start - ex.y0)/ex.dy)
    j_f = int((y_stop - ex.y0)/ex.dy)
    return grid_x[i_0:i_f], grid_y[j_0:j_f]
    
#--------------------------

def plot_bfs_compare_N(Ns, N_max):
    err_xs, err_ys = cnvg.bfs_compare_N(example, Ns, N_max)
    title_xs = "Error to $N^{*}=$%d in reattatchment point"%N_max
    title_ys = "Error to $N^{*}=$%d in detachment point"%N_max
    ax_labels_detach= ["N", "$|y_{N^{*}} - y_{N}|$"]
    ax_labels_reattach= ["N", "$|x_{N^{*}} - x_{N}|$"]
    
    n_feats = 3
    
    labels_stream = np.arange(1, n_feats+1)
    
    graphics.plot_log_multi(err_xs[:n_feats], Ns, title_xs, labels_stream, ax_labels_reattach)
    graphics.plot_log_multi(err_ys[:n_feats], Ns, title_ys, labels_stream, ax_labels_detach)

def plot_tri_compare_N(Ns, N_max):
    err_extrs, err_extrs_y, err_left, err_right = cnvg.tri_compare_N(example, Ns, N_max)
    
    title_extrs_y = "Error to $N^{*}=$%d in $y$ of vortex center"%N_max
    title_extrs = "Error to $N^{*}=$%d in $\psi$ of vortex center"%N_max
    ax_labels_stream = ["N", "$|\psi_{N^{*}} - \psi_{N}|$"]
    ax_labels_y = ["N", "$|y_{N^{*}} - y_{N}|$"]
    
    n_feats = 3
    
    labels_stream_extrs = np.arange(1, n_feats+1)
    
    graphics.plot_log_multi(err_extrs_y[:n_feats], Ns, title_extrs_y, labels_stream_extrs, ax_labels_y)
    graphics.plot_log_multi(err_extrs[:n_feats], Ns, title_extrs, labels_stream_extrs, ax_labels_stream)

    
    title_left = "Error to $N^{*}=$%d in $y$ of stream-saddle on left boundary"%N_max
    title_right = "Error to $N^{*}=$%d in $y$ of stream-saddle on right boundary"%N_max
    ax_labels_saddle = ["N", "$|y_{N^{*}} - y_{N}|$"]
    
    labels_stream_left = np.arange(1, n_feats+1)
    labels_stream_right = np.arange(1, n_feats+1)
    
    graphics.plot_log_multi(err_left[:n_feats], Ns, title_left, labels_stream_left, ax_labels_saddle)
    graphics.plot_log_multi(err_right[:n_feats], Ns, title_right, labels_stream_right, ax_labels_saddle)