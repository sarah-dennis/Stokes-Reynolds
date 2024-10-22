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
import stokes_pressure as pressure

#------------------------------------------------------------------------------
from stokes_solver import run_spLU

class Stokes_Solver:
    def __init__(self, example, max_iters=50000):
        self.example = example
        self.max_iters = max_iters
        self.write_mod = 200
        self.error_mod = 200
        self.err_tol = 1e-8
#------------------------------------------------------------------------------
    def new_run(self, N, iters):
        ex = self.example(N)
        
        u_init = np.zeros(ex.Nx * ex.Ny)
        v_init = np.zeros(ex.Nx * ex.Ny)
        psi_init = np.zeros(ex.Nx * ex.Ny)
        past_iters = 0
    
        u, v, psi = run_spLU(ex, u_init, v_init, psi_init, iters, past_iters, self.error_mod, self.write_mod, self.err_tol)
    
        rw.write_solution(ex, u, v, psi, iters)
    
    def new_run_many(self, N_0, dN, many):
        
        self.new_run(N_0, self.max_iters)
        N_load = N_0
        for k in range (1, many): 
            N = N_load*dN
            self.load_scale(N_load, N)
            self.load_run(N, self.max_iters)
            N_load = N

    
    def load_run_new_many(self, N_0, dN, many):
        N_load = N_0
        self.load_run(N_load, self.max_iters)
        for k in range (many): 
            N = N_load*dN
            self.load_scale(N_load, N)
            self.load_run(N, self.max_iters)
            N_load = N

    def load_run_many(self, N_0, dN, many):
        N = N_0
        for k in range (many): 
            self.load_run(N, self.max_iters)
            N *= dN 
                                                                                                                                                                                                                                                                       
    def load_run(self, N, iters):                                
        ex = self.example(N)
        u, v, psi, past_iters = rw.read_solution(ex.filestr+".csv", ex.Nx*ex.Ny)
        
        u, v, psi = run_spLU(ex, u, v, psi, iters, past_iters, self.error_mod, self.write_mod, self.err_tol)
        
        rw.write_solution(ex, u, v, psi, iters+past_iters)
    

    def load_scale(self, N_load, N_scale):
        ex_load = self.example(N_load)
        ex_scale = self.example(N_scale)
        
        points_load = (ex_load.ys, ex_load.xs)
        
        u_load, v_load, psi_load, past_iters = rw.read_solution(ex_load.filestr+".csv", ex_load.Ny*ex_load.Nx)
        u_load_2D = u_load.reshape((ex_load.Ny,ex_load.Nx))
        v_load_2D = v_load.reshape((ex_load.Ny,ex_load.Nx))
        psi_load_2D = psi_load.reshape((ex_load.Ny,ex_load.Nx))
    
    
        points_scale = np.meshgrid(ex_scale.ys, ex_scale.xs)
        
        u_scaled_2D = interpn(points_load, u_load_2D, tuple(points_scale))
        v_scaled_2D = interpn(points_load, v_load_2D, tuple(points_scale))
        psi_scaled_2D = interpn(points_load, psi_load_2D, tuple(points_scale))
    
        u_scaled = u_scaled_2D.ravel(order='F')
        v_scaled = v_scaled_2D.ravel(order='F')
        psi_scaled = psi_scaled_2D.ravel(order='F')
    
        rw.write_solution(ex_scale, u_scaled, v_scaled, psi_scaled, 0)
    
#------------------------------------------------------------------------------
# Error
#------------------------------------------------------------------------------
    def compare(self, N_min, Ns, N_max):
        l1_errs, l2_errs, inf_errs, cnvg_rates, ex_min = cnvg.compare_Ns(self.example, N_min, Ns, N_max)
        title = "Error in $\psi$ to $N^{*}=%d$ \n %s"%(N_max, ex_min.spacestr)
        ax_labels = ["N", "$||\psi _{N^{*}} - \psi_{N}||_p$"]
        leg_labels = ['$L^1$', '$L^2$','$L^\infty$', 'res_err']
        linthresh = 1e-8
        O1 = 1
        O2 = 1
        graphics.plot_log_multi([l1_errs, l2_errs, inf_errs], [N_min]+Ns, title, leg_labels, ax_labels, linthresh, O1, O2)
        print("cnvg rates")
        print("l1: " + np.array2string(cnvg_rates[0]))
        print("l2: " + np.array2string(cnvg_rates[1]))
        print("linf" + np.array2string(cnvg_rates[2]))
    
#------------------------------------------------------------------------------
# PLOTTING 
#------------------------------------------------------------------------------
    def load_plot(self, N):
        ex = self.example(N)
        u, v, psi, past_iters = rw.read_solution(ex.filestr+".csv", ex.Nx * ex.Ny)
    
    # Grid domain
        xs = ex.xs
        ys = ex.ys
        graphics.plot_contour_mesh(ex.space, xs, ys, 'space',['space', 'x', 'y'],log_cmap=False)
    
    
    
    
    # zoom domain for vorticity & pressure
        x_start = 0.9
        x_stop= 1.1
        y_start = 0.9
        y_stop = 1.1
        xs_zoom, ys_zoom = grid_zoom_1D(xs, ys, ex, x_start, x_stop, y_start, y_stop)

    # Pressure plot: 
        p = pressure.pressure(ex, u, v)
        p_2D = p.reshape((ex.Ny,ex.Nx))
        dp, res = pressure.resistance(ex, p) 
        dp_str = ', $\Delta P =%.2f$'%dp
    
        p_max = np.max(p)
        p_min = np.min(p)
    
        ax_labels_p = ['$p(x,y)$', '$x$', '$y$']
        title_p = 'Pressure $p(x,y)$ \n' + ex.spacestr + dp_str
    
        p_ma = np.ma.masked_where(ex.space==-1, p_2D)
    
        graphics.plot_contour_mesh(p_ma, xs, ys, title_p, ax_labels_p, log_cmap=False , n_contours=40, vmax=p_max, vmin=p_min)
    
        p_zoom = grid_zoom_2D(p_ma, ex, x_start, x_stop, y_start, y_stop)     
        graphics.plot_contour_mesh(p_zoom, xs_zoom, ys_zoom, title_p, ax_labels_p, log_cmap=False, n_contours=20, vmax=p_max, vmin=p_min)
    
    # zoom domain for velocity & stream
        x_start = 1
        x_stop= 1.5
        y_start = 0
        y_stop = 0.5
        xs_zoom, ys_zoom = grid_zoom_1D(xs, ys, ex, x_start, x_stop, y_start, y_stop)
    
    #  Velocity plot: 
    
        ax_labels = ['$|(u,v)|_2$','$x$', '$y$']
        title = 'Velocity $(u,v)$ \n' + ex.spacestr + dp_str
        ax_labels = ['$|(u,v)|_2$','$x$', '$y$']
        
        u_2D = u.reshape((ex.Ny,ex.Nx))
        v_2D = v.reshape((ex.Ny,ex.Nx))
    
        uv_mag = np.sqrt(u_2D**2 + v_2D**2)
        uv_mag_max = np.max(uv_mag)
        
        u_2D_ma = np.ma.masked_where(ex.space==-1,u_2D)
        v_2D_ma = np.ma.masked_where(ex.space==-1,v_2D)
        graphics.plot_stream_heat(u_2D_ma, v_2D_ma, xs, ys, uv_mag, title, ax_labels, log_cmap=False, vmin=0, vmax=uv_mag_max) 
        
        u_2D_zoom = grid_zoom_2D(u_2D_ma, ex, x_start, x_stop, y_start, y_stop)
        v_2D_zoom = grid_zoom_2D(v_2D_ma, ex, x_start, x_stop, y_start, y_stop)
        uv_mag_zoom = grid_zoom_2D(uv_mag, ex, x_start, x_stop, y_start, y_stop)
        graphics.plot_stream_heat(u_2D_zoom, v_2D_zoom, xs_zoom, ys_zoom, uv_mag_zoom, title, ax_labels, log_cmap=False, vmin=0, vmax=uv_mag_max)
    
    # Stream plot:
    
        # ax_labels = ['$\psi(x,y)$ : $u = \psi_y$, $v = -\psi_x$', '$x$', '$y$']
        # title = 'Stream $\psi(x,y)$ \n' + ex.spacestr + dp_str
        # stream_2D = psi.reshape((ex.Ny,ex.Nx))
        # stream_2D_ma = np.ma.masked_where(ex.space==-1, stream_2D)
        # graphics.plot_contour_mesh(stream_2D_ma, xs, ys, title, ax_labels, log_cmap=True, n_contours=20, vmin=None, vmax=ex.flux)
    
        # stream_2D_zoom = grid_zoom_2D(stream_2D_ma, ex, x_start, x_stop, y_start, y_stop)
        # graphics.plot_contour_mesh(stream_2D_zoom, xs_zoom, ys_zoom, title, ax_labels, True, n_contours=20, vmin=None, vmax=ex.flux)
    

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
    
