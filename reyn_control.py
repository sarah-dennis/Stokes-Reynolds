# -*- coding: utf-8 -*-
"""
Created on Tue May 21 16:05:05 2024

@author: sarah
"""
import numpy as np
import graphics
import readwrite as rw
import convergence as cnvg

import reyn_velocity as rv
import reyn_pressure as rp 
import reyn_perturbed as rpert



lenx = 2
leny = 1
x_start = 0
y_start = 0
x_stop= x_start + lenx
y_stop = y_start + leny


log_linthresh=1e-8  

class Reynolds_Solver: 
    def __init__(self, Example, U, dP, args=None):
        self.Example = Example #initialize ex = Example(U, dP, args) in the solver
        self.args=args
        self.U = U
        self.dP = dP

        # colorbar min max

        self.vel_max = 5
        self.p_min=-3
        self.p_max=3

    

    
    def pwl_solve(self, N, write=False, plot=True,zoom=False):
        ex = self.Example(self.U, self.dP, N, self.args)

        pressure = rp.PwlGMRes_ReynPressure(ex)
        velocity = rv.ReynVelocity(ex, pressure.ps_1D)
        solver_title = "Reynolds Piecewise Linear"
        if plot:
            self.p_plot(ex, pressure, velocity.flux, solver_title, zoom)
            self.v_plot(ex, velocity, pressure.dP, solver_title, zoom)
        
        if write:
            nm = ex.Nx * ex.Ny
            u = velocity.vx.reshape(nm)
            v = velocity.vy.reshape(nm)
            p = pressure.ps_2D.reshape(nm)
            rw.write_reyn(ex, u, v, p)
        
        return pressure, velocity
    
    def fd_solve(self, N, write=False, plot=True,zoom=False):
        ex = self.Example(self.U, self.dP, N, self.args)
        pressure = rp.FinDiff_ReynPressure(ex)
        
        velocity = rv.ReynVelocity(ex, pressure.ps_1D)
        solver_title = "Reynolds Finite Difference"
        if plot:
            self.p_plot(ex, pressure, velocity.flux, solver_title, zoom)
            self.v_plot(ex, velocity, pressure.dP, solver_title, zoom)

        if write:
            nm = ex.Nx * ex.Ny
            u = velocity.vx.reshape(nm)
            v = velocity.vy.reshape(nm)
            p= pressure.ps_2D.reshape(nm)
            rw.write_reyn(ex, u, v, p)
        return pressure, velocity
    
    def fd_adj_solve(self, N, write=False, plot=True, zoom=False):
        ex = self.Example(self.U, self.dP, N, self.args)
        
        adj_pressure = rp.Adjusted_ReynPressure(ex)
        adj_velocity = rv.Adjusted_ReynVelocity(ex, adj_pressure.ps_2D)
        solver_title = "Adjusted Reynolds"
        if plot:
            self.p_plot(ex, adj_pressure , adj_velocity.flux, solver_title, zoom)
            self.v_plot(ex, adj_velocity, adj_pressure.dP, solver_title,zoom)
        
        if write:
            nm = ex.Nx * ex.Ny
            u = adj_velocity.vx.reshape(nm)
            v = adj_velocity.vy.reshape(nm)
            p = adj_pressure.ps_2D.reshape(nm)
            rw.write_reyn(ex, u, v, p)
        
        return adj_pressure, adj_velocity
    
    

    def fd_pert_solve(self, N, order, write=False, plot=True, zoom=False, get_dPs = False):
        ex = self.Example(self.U, self.dP, N, self.args)

        
        reyn_pressure = rp.FinDiff_ReynPressure(ex)
        reyn_velocity = rv.ReynVelocity(ex, reyn_pressure.ps_1D)                   
        
        pert = rpert.PerturbedReynSol(ex, order, reyn_pressure, reyn_velocity)
        solver_title = "Reynolds"
        
        if plot:
            # self.p_plot(ex, reyn_pressure, reyn_velocity.flux, solver_title, zoom)
            # self.v_plot(ex, reyn_velocity, reyn_pressure.dP, solver_title, zoom)
            if order > 1:
                solver_title2 = solver_title + " $O(\delta^2)$ perturbed"
                self.p_plot(ex, pert.pert2_pressure, pert.pert2_velocity.flux, solver_title2, zoom)
                self.v_plot(ex, pert.pert2_velocity, pert.pert2_pressure.dP, solver_title2, zoom)
            if order > 3:
                solver_title4 = solver_title + " $O(\delta^4)$ perturbed"
                self.p_plot(ex, pert.pert4_pressure, pert.pert4_velocity.flux, solver_title4, zoom)
                self.v_plot(ex, pert.pert4_velocity, pert.pert4_pressure.dP, solver_title4, zoom)
        
        if write:
            if order > 1 and order < 3:
                nm = ex.Nx * ex.Ny
                u = pert.pert2_velocity.vx.reshape(nm)
                v = pert.pert2_velocity.vy.reshape(nm)
                p = pert.pert2_pressure.ps_2D.reshape(nm)
                rw.write_reyn(ex, u, v, p)
            if order > 1 and order < 5:
                nm = ex.Nx * ex.Ny
                u = pert.pert4_velocity.vx.reshape(nm)
                v = pert.pert4_velocity.vy.reshape(nm)
                p = pert.pert4_pressure.ps_2D.reshape(nm)
                rw.write_reyn(ex, u, v, p)
            
        if get_dPs:
            return pert
        else:
            if order <3:
                return pert.pert2_pressure, pert.pert2_velocity
            else:
                return pert.pert4_pressure, pert.pert4_velocity
    
    
    def p_plot(self, ex, pressure, flux, solver_title, zoom):
        paramstr = "$Re=0$, $Q=%.2f$, $U=%.1f$, $\Delta P=%.2f$"%(flux, self.U, pressure.dP)
        p_title = solver_title +'\n' + paramstr
        p_labels = ["$p(x)$", "$x$","$y$"]
        graphics.plot_contour_mesh(pressure.ps_2D, ex.xs, ex.ys, p_title, p_labels, vmin=self.p_min, vmax=self.p_max, log_cmap=False)
    
        if zoom:
            xs_zoom, ys_zoom = graphics.grid_zoom_1D(ex.xs, ex.ys, ex, x_start, x_stop, y_start, y_stop)
            p_zoom = graphics.grid_zoom_2D(pressure.ps_2D, ex, x_start, x_stop, y_start, y_stop)
            graphics.plot_contour_mesh(p_zoom, xs_zoom, ys_zoom, p_title, p_labels, vmin=self.p_min, vmax=self.p_max, log_cmap=False)
    
    
    def v_plot(self, ex, velocity, dP, solver_title, zoom):
        paramstr = "$Re=0$, $Q=%.2f$, $U=%.2f$, $\Delta P=%.2f$"%(velocity.flux, self.U, dP)
        v_title = solver_title + '\n' + paramstr
        v_ax_labels =  ['$|(u,v)|_2$','$x$', '$y$'] 
        uv_mag = np.sqrt(velocity.vx**2 + velocity.vy**2)
        
        graphics.plot_stream_heat(velocity.vx, velocity.vy, ex.xs, ex.ys, uv_mag, v_title, v_ax_labels, vmin=0, vmax=self.vel_max, log_cmap=False)

        # graphics.plot_contour_mesh(uv_mag, ex.xs, ex.ys, v_title, v_ax_labels, vmin=0, vmax=self.vel_max, log_cmap=False)
        # graphics.plot_quiver(velocity.vx, velocity.vy, ex.xs, ex.ys, uv_mag, v_title, v_ax_labels, vmin=0, vmax=self.vel_max)

        # graphics.plot_contour_mesh(velocity.vx, ex.xs, ex.ys, 'u', ['$u$', '$x$', '$y$'], -3, 3)
        # graphics.plot_contour_mesh(velocity.vy, ex.xs, ex.ys, 'v', ['$v$', '$x$', '$y$'], -3, 3)
        
        graphics.plot_contour_mesh(velocity.inc, ex.xs, ex.ys, 'incompressibility', ['$u_x+v_y$', '$x$', '$y$'], -1/1, 1/1)
        
        graphics.plot_2D(velocity.qs, ex.xs, 'flux $\mathcal{Q} = q(x) =\int_0^{h(x)} u(x,y) dy$', ['$x$', '$q(x)=\mathcal{Q}$'])

        i_test = 3*ex.Nx//4
        graphics.plot_2D(ex.ys,velocity.vx[:,i_test],  f'$u({ex.xs[i_test]:.2f}),y)$',[f'$u({ex.xs[i_test]:.2f},y)$','$y$'])
        graphics.plot_2D(ex.ys,velocity.vy[:,i_test],  f'$v({ex.xs[i_test]:.2f}),y)$',[f'$v({ex.xs[i_test]:.2f},y)$','$y$'])
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
        pressure = rp.Pressure(ex, ps_2D=p)
        velocity = rv.Velocity(ex, u, v)
        self.p_plot(ex, pressure, velocity.flux, zoom)
        self.v_plot(ex, velocity, zoom)

    def convg_pwl_fd(self, Ns):
        l1_errs, l2_errs, inf_errs, cnvg_rates = cnvg.reyn_cnvg_pwl_fd(self, Ns)
        title ='' #'Grid Convergence for Reynolds pressure \n Finite-Difference to Piecewise-Linear-Analytic'
        ax_labels=['$N$',  '$||p_{FD} - p_{PLA}||$']
        fun_labels=['$L_1$', '$L_2$', '$L_\infty$']
        O1 = 1
        O2 = 1
        graphics.plot_log_multi([l1_errs, l2_errs,inf_errs],Ns,title,fun_labels,ax_labels,log_linthresh, O1, O2)
   