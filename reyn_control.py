# -*- coding: utf-8 -*-
"""
Created on Tue May 21 16:05:05 2024

@author: sarah
"""
import numpy as np
import graphics
import readwrite as rw
import convergence as cnvg
import reyn_pressures_pwlinear as pwl


import reyn_pressures_finDiff as fd

from reyn_velocity import Velocity, Adj_Velocity
from reyn_pressure import Pressure, Adj_Pressure
from numpy.linalg import solve as np_fd_solve
from reyn_heights import PWL_Height
from scipy.sparse.linalg import gmres as sp_pwl_gmres


lenx = 1
leny = 1
x_start = 1
y_start =0.5
x_stop= x_start + lenx
y_stop = y_start + leny


log_linthresh=1e-8  

class Reynolds_Solver: 
    def __init__(self, Example, U, dP, args=None):
        self.Example = Example #initialize ex = Example(U, dP, args) in the solver
        self.args=args
        self.U = U
        self.dP = dP

        # colorbar min maxs
        self.vel_max =70
        self.p_min=-10
        self.p_max=10
    

    
    def pwl_solve(self, N, plot=True,zoom=False):
        ex = self.Example(self.U, self.dP, N, self.args)
        
        if not isinstance(ex, PWL_Height):
            raise TypeError('Example is not piecewise linear')
        
        rhs = pwl.make_rhs(ex)
        linOp = pwl.pwlLinOp(ex)
        coefs, exit_code = sp_pwl_gmres(linOp, rhs, tol=1e-10)
         
        if exit_code != 0:
            raise Exception('gmres did not converge')
        ps, flux = pwl.make_ps(ex, coefs)
        pressure = Pressure(ex, ps)
        velocity = Velocity(ex, ps)
        
        
        if plot:
            self.p_plot(ex, pressure, flux,zoom)
            
            self.v_plot(ex, velocity, zoom)
        return pressure, velocity
    
    def fd_solve(self, N, plot=True,zoom=False):
        ex = self.Example(self.U, self.dP, N, self.args)
        rhs = fd.make_rhs(ex)
        mat = fd.make_mat(ex)
        ps = np_fd_solve(mat, rhs)
        pressure  = Pressure(ex, ps)
        velocity = Velocity(ex, ps)
        if plot:
            self.p_plot(ex, pressure, velocity.flux, zoom)
            
            self.v_plot(ex, velocity, zoom)
        return pressure, velocity
    
    def fd_adj_solve(self, N, plot=True, zoom=False):
        ex = self.Example(self.U, self.dP, N, self.args)
        
        rhs = fd.make_rhs(ex)
        mat = fd.make_mat(ex)
        p_reyn = np_fd_solve(mat, rhs)
        adj_pressure = Adj_Pressure(ex, p_reyn)
        adj_velocity = Adj_Velocity(ex, adj_pressure.ps_2D)
        
        if plot:
            self.p_plot(ex, adj_pressure, adj_velocity.flux, zoom)
            self.v_plot(ex, adj_velocity, zoom)
        
        nm = ex.Nx * ex.Ny
        u = adj_velocity.vx.reshape(nm)
        v = adj_velocity.vy.reshape(nm)
        p_adj = adj_pressure.ps_2D.reshape(nm)
        rw.write_reyn(ex, u, v, p_adj)
        
        return adj_pressure, adj_velocity

    def p_plot(self, ex, pressure, flux, zoom):
        paramstr = "$Re=0$, $Q=%.2f$, $U=%.1f$, $\Delta P=%.2f$"%(flux, self.U, self.dP)
        p_title = "" + paramstr
        p_labels = ["$p(x)$", "$x$","$y$"]
        graphics.plot_contour_mesh(pressure.ps_2D, ex.xs, ex.ys, p_title, p_labels, vmin=self.p_min, vmax=self.p_max, log_cmap=False)
    
        if zoom:
            xs_zoom, ys_zoom = graphics.grid_zoom_1D(ex.xs, ex.ys, ex, x_start, x_stop, y_start, y_stop)
            p_zoom = graphics.grid_zoom_2D(pressure.ps_2D, ex, x_start, x_stop, y_start, y_stop)
            graphics.plot_contour_mesh(p_zoom, xs_zoom, ys_zoom, p_title, p_labels, vmin=self.p_min, vmax=self.p_max, log_cmap=False)
    
    
    def v_plot(self, ex, velocity, zoom):
        paramstr = "$Re=0$, $Q=%.2f$, $U=%.1f$, $\Delta P=%.2f$"%(velocity.flux, self.U, self.dP)
        v_title = "" + paramstr
        v_ax_labels =  ['$|(u,v)|_2$','$x$', '$y$'] 
        uv_mag = np.sqrt(velocity.vx**2 + velocity.vy**2)
        graphics.plot_stream_heat(velocity.vx, velocity.vy, ex.xs, ex.ys, uv_mag, v_title, v_ax_labels, vmin=0, vmax=self.vel_max, log_cmap=False)

        if zoom:

            xs_zoom, ys_zoom = graphics.grid_zoom_1D(ex.xs, ex.ys, ex, x_start, x_stop, y_start, y_stop)
            u_2D_zoom = graphics.grid_zoom_2D(velocity.vx, ex, x_start, x_stop, y_start, y_stop)
            v_2D_zoom = graphics.grid_zoom_2D(velocity.vy, ex, x_start, x_stop, y_start, y_stop)
            uv_mag_zoom = graphics.grid_zoom_2D(uv_mag, ex, x_start, x_stop, y_start, y_stop)
            graphics.plot_stream_heat(u_2D_zoom, v_2D_zoom, xs_zoom, ys_zoom, uv_mag_zoom, v_title, v_ax_labels, vmin=0, vmax=self.vel_max, log_cmap=False)

    def load_plot(self, N, zoom):
        ex = self.Example(self.U, self.dP, N, self.args)
        u, v, p = rw.read_reyn(ex.filestr+'.csv', ex.Nx, ex.Ny)
        u = u.reshape((ex.Ny, ex.Nx))
        v = v.reshape((ex.Ny, ex.Nx))
        p = p.reshape((ex.Ny, ex.Nx))
        pressure = Adj_Pressure(ex, p)
        velocity = Adj_Velocity(ex, p)
        self.p_plot(ex, pressure, velocity.flux, zoom)
        self.v_plot(ex, velocity, zoom)

    def compare_pwl_fd(self, Ns):
        
        l1_errs, l2_errs, inf_errs, cnvg_rates = cnvg.reyn_cnvg_pwl_fd(self, Ns)
        
        
        title ='' #'Grid Convergence for Reynolds pressure \n Finite-Difference to Piecewise-Linear-Analytic'
        ax_labels=['$N$',  '$||p_{FD} - p_{PLA}||$']
        fun_labels=['$L_1$', '$L_2$', '$L_\infty$']
        O1 = 1
        O2 = 1
        graphics.plot_log_multi([l1_errs, l2_errs,inf_errs],Ns,title,fun_labels,ax_labels,log_linthresh, O1, O2)
    




