#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 22 14:59:46 2024

@author: sarahdennis
"""
import numpy as np
from scipy.interpolate import interpn

import graphics
import readwrite as rw
import convergence as cnvg
import stokes_pressure as pressure

#------------------------------------------------------------------------------
from stokes_solver import run_spLU
#---------------------------------------------------------------------------
# graphics args

#-> zoom plot
lenx = 0.55
leny = 0.55
x_start = 2
y_start = 1.5
x_stop= x_start + lenx
y_stop = y_start + leny

#-> log plots
log_linthresh=1e-5  
log_cmap_on = False

#---------------------------------------------------------------------------
class Stokes_Solver:
    def __init__(self, Example, args, U, Q, Re, max_iters=50000):
        # domain 
        self.Example = Example
        self.args = args
        self.U = U
        self.Q = Q
        self.Re = Re
        
        # iterative solution args
        self.max_iters = max_iters
        self.write_mod = 500
        self.error_mod = 500
        self.err_tol = 1e-8
        
        # plotting thresholds

        self.vel_max = 5
        self.p_min=-0
        self.p_max=40

        
#------------------------------------------------------------------------------
    def new_run(self, N):
        ex = self.Example(self.args, self.U, self.Q, self.Re, N)
        
        u_init = np.zeros(ex.Nx * ex.Ny)
        v_init = np.zeros(ex.Nx * ex.Ny)
        psi_init = np.zeros(ex.Nx * ex.Ny)
        past_iters = 0
    
        u, v, psi = run_spLU(ex, u_init, v_init, psi_init, self.max_iters, past_iters, self.error_mod, self.write_mod, self.err_tol)
    
        rw.write_stokes(ex, u, v, psi, self.max_iters)
    
    def new_run_many(self, N_0, dN, many):
        self.new_run(N_0)
        N_load = N_0
        for k in range (1, many): 
            N = N_load*dN
            self.load_scale(N_load, N)
            self.load_run(N)
            N_load = N

    def load_run_new_many(self, N_0, dN, many):
        N_load = N_0
        self.load_run(N_load)
        for k in range (many): 
            N = N_load*dN
            self.load_scale(N_load, N)
            self.load_run(N)
            N_load = N

    def load_run_many(self, N_0, dN, many):
        N = N_0
        for k in range (many): 
            self.load_run(N)
            N *= dN 
                                                                                                                                                                                                                                                                       
    def load_run(self, N):                                
        ex = self.Example(self.args, self.U, self.Q, self.Re, N)
        u, v, psi, past_iters = rw.read_stokes(ex.filestr+".csv", ex.Nx*ex.Ny)
        u, v, psi = run_spLU(ex, u, v, psi, self.max_iters, past_iters, self.error_mod, self.write_mod, self.err_tol)
        rw.write_stokes(ex, u, v, psi, self.max_iters+past_iters)

    def load_scale(self, N_load, N_scale):
        ex_load = self.Example(self.args, self.U, self.Q, self.Re, N_load)
        ex_scale = self.Example(self.args, self.U, self.Q, self.Re, N_scale)
        
        points_load = (ex_load.ys, ex_load.xs)
        
        u_load, v_load, psi_load, past_iters = rw.read_stokes(ex_load.filestr+".csv", ex_load.Ny*ex_load.Nx)
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
    
        rw.write_stokes(ex_scale, u_scaled, v_scaled, psi_scaled, 0)
    
    def load(self,N):
        ex = self.Example(self.args, self.U, self.Q, self.Re, N)
        u, v, psi, past_iters = rw.read_stokes(ex.filestr+".csv", ex.Nx*ex.Ny)
        p = pressure.pressure(ex, u, v)
        p = p.reshape((ex.Ny,ex.Nx))
        u = u.reshape((ex.Ny,ex.Nx))
        v = v.reshape((ex.Ny,ex.Nx))
        
        p = np.flip(p, axis=0)        
        u = np.flip(u, axis=0)
        v = -np.flip(v, axis=0) 
        
        return p, u, v
        
#------------------------------------------------------------------------------    
    def get_dP(self,N):
        ex = self.Example(self.args, self.U, self.Q, self.Re, N)
        u, v, psi, past_iters = rw.read_stokes(ex.filestr+".csv", ex.Nx * ex.Ny)
        p = pressure.pressure(ex, u, v)
        dp, res = pressure.resistance(ex, p) 
        return dp
    
    def get_attachments(self,N):
        ex=self.Example(self.args, self.U, self.Q, self.Re, N)
        u, v, psi, past_iters = rw.read_stokes(ex.filestr+".csv", ex.Nx * ex.Ny)
        
        y_xr = 0 
        x_yr = 1 #=l
        
        i_yr =int( x_yr/ex.dx +1)
        
        j_xr = int(y_xr+1)

        xr = 0
        yr = 0
        
        for i in range(ex.Nx-1):
            k = j_xr*ex.Nx + i
            k2 = j_xr*ex.Nx + i+1    
            if np.sign(psi[k])!= np.sign(psi[k2]):
                xr = ex.xs[i]-1
                
        for j in range(ex.Ny-1):
            k = int(j*ex.Nx + i_yr)
            k2 = int((j+1)*ex.Nx + i_yr)
            if np.sign(psi[k])!= np.sign(psi[k2]):
                yr = ex.ys[j]
        return xr,yr
    
#------------------------------------------------------------------------------
# Error
#------------------------------------------------------------------------------
    def compare(self,N_min, Ns, N_max,p_err=True):
        
        l1_errs, l2_errs, inf_errs, cnvg_rates, ex_min = cnvg.stokes_cnvg_self(self.Example, N_min, Ns, N_max,p_err)
        title = ''#"Iterative Grid Error in Stream $\psi$ at $N_{max}=%d$ \n %s"%(N_max, ex_min.spacestr)
        ax_labels = ["$N$", "$||\psi _{N^{*}} - \psi_{N}||_p$"]
        leg_labels = ['$L^1$', '$L^2$','$L^\infty$']
        
        O1 = 1
        O2 = 1
        graphics.plot_log_multi([l1_errs, l2_errs, inf_errs], [N_min]+Ns, title, leg_labels, ax_labels, log_linthresh, O1, O2)

#------------------------------------------------------------------------------
# PLOTTING 
#------------------------------------------------------------------------------
    def load_plot(self, N,zoom=False):
        ex = self.Example(self.args, self.U, self.Q, self.Re, N)
        u, v, psi, past_iters = rw.read_stokes(ex.filestr+".csv", ex.Nx * ex.Ny)

    
    # Grid domain
        xs = ex.xs
        ys = ex.ys
        
    # Space plot (xi,yj: boundary/interior/exterior)
        # graphics.plot_contour_mesh(ex.space, xs, ys, 'space',['space', 'x', 'y'],log_cmap=False)
    

    # Zoom domain
        if zoom:
            xs_zoom, ys_zoom = graphics.grid_zoom_1D(xs, ys, ex, x_start, x_stop, y_start, y_stop)

    # Pressure plot: 
        p = pressure.pressure(ex, u, v)
        dp, res = pressure.resistance(ex, p) 
        
        p_2D = p.reshape((ex.Ny,ex.Nx))
        dp_str = ', $\Delta P =%.2f$'%(-dp)

    
        ax_labels_p = ['$p(x,y)$', '$x$', '$y$']
        title_p = 'Stokes\n' + ex.spacestr + dp_str
    
        p_ma = np.ma.masked_where(ex.space==-1, p_2D)
        p_ma = np.flip(p_ma, axis=0)
        graphics.plot_contour_mesh(p_ma, xs, ys, title_p, ax_labels_p,  vmin=self.p_min, vmax=self.p_max, log_cmap=log_cmap_on, n_contours=40)#, y_lim = np.min(ex.y_peaks))
    
        if zoom:
            p_zoom = graphics.grid_zoom_2D(p_ma, ex, x_start, x_stop, y_start, y_stop)  
            graphics.plot_contour_mesh(p_zoom, xs_zoom, ys_zoom, title_p, ax_labels_p, vmin=self.p_min, vmax=self.p_max, log_cmap=log_cmap_on, n_contours=20)#, y_lim = min(ex.y_peaks))
    
    #  Velocity plot: 
    
        ax_labels = ['$|(u,v)|_2$','$x$', '$y$']
        
        title = 'Stokes\n' + ex.spacestr + dp_str
        ax_labels = ['$|(u,v)|_2$','$x$', '$y$']
        
        u_2D = u.reshape((ex.Ny,ex.Nx))

        v_2D = v.reshape((ex.Ny,ex.Nx)) 

        uv_mag = np.sqrt(u_2D**2 + v_2D**2)
        
        
        uv_mag = np.flip(uv_mag, axis=0)
        u_2D = np.flip(u_2D, axis=0)
        v_2D = -np.flip(v_2D, axis=0)
        
        graphics.plot_stream_heat(u_2D, v_2D, xs, ys, uv_mag, title, ax_labels,  vmin=0, vmax=self.vel_max, log_cmap=False) 
        
        if zoom:
            u_2D_zoom = graphics.grid_zoom_2D(u_2D, ex, x_start, x_stop, y_start, y_stop)
            v_2D_zoom = graphics.grid_zoom_2D(v_2D, ex, x_start, x_stop, y_start, y_stop)
            uv_mag_zoom = graphics.grid_zoom_2D(uv_mag, ex, x_start, x_stop, y_start, y_stop)
            graphics.plot_stream_heat(u_2D_zoom, v_2D_zoom, xs_zoom, ys_zoom, uv_mag_zoom, title, ax_labels, vmin=0, vmax=self.vel_max, log_cmap=False)
        
    # Stream plot:
    
        # ax_labels = ['$\psi(x,y)$ : $u = \psi_y$, $v = -\psi_x$', '$x$', '$y$']
        # title = 'Stream $\psi(x,y)$ \n' + ex.spacestr + dp_str
        # stream_2D = psi.reshape((ex.Ny,ex.Nx))
        # stream_2D_ma = np.ma.masked_where(ex.space==-1, stream_2D)
        # graphics.plot_contour_mesh(stream_2D_ma, xs, ys, title, ax_labels, log_cmap=True, n_contours=20, vmin=None, vmax=ex.flux)
        # if zoom:
            # stream_2D_zoom = grid_zoom_2D(stream_2D_ma, ex, x_start, x_stop, y_start, y_stop)
            # graphics.plot_contour_mesh(stream_2D_zoom, xs_zoom, ys_zoom, title, ax_labels, True, n_contours=20, vmin=None, vmax=ex.flux)
    
        return dp
