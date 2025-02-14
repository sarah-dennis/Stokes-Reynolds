# -*- coding: utf-8 -*-
"""
Created on Tue May 21 16:05:05 2024

@author: sarah
"""
import numpy as np
import graphics
import reyn_pressures_pwlinear as pwl
from scipy.sparse.linalg import gmres as sp_gmres

import reyn_pressures_finDiff as fd
from numpy.linalg import solve as np_solve

from reyn_velocity import Velocity, Adj_Velocity
from reyn_pressure import Pressure, Adj_Pressure


class Reynolds_Solver: 
    def __init__(self, Example, U, dP, args=None):
        self.Example = Example #initialize ex = Example(U, dP, args) in the solver
        self.args=args
        self.U = U
        self.dP = dP

        # colorbar min maxs
        self.vel_max = 4
        self.p_min=-70
        self.p_max=70
    
    def est_flux(self,height,ps):
        i=height.Nx//2
        j=height.Ny-5
        dp = (ps[j][i+1]-ps[j][i-1])/(2*height.dx)
        hi=height.hs[i]
        flux = -(hi**3 * dp)/(12*height.visc) + height.U*hi/2
        return flux
        
    
    def pwl_solve(self, N, plot=True):
        ex = self.Example(self.U, self.dP, N, self.args)
        rhs = pwl.make_rhs(ex)
        linOp = pwl.pwlLinOp(ex)
        coefs, exit_code = sp_gmres(linOp, rhs, tol=1e-10)
         
        if exit_code != 0:
            print('gmres did not converge')
        ps, flux = pwl.make_ps(ex, coefs)
        pressure = Pressure(ex, ps)
        
        if plot:
            self.p_plot(ex, pressure, flux)
            velocity = Velocity(ex, ps)
            self.v_plot(ex, velocity, flux)
        return pressure, flux
    
    def fd_solve(self, N, plot=True):
        ex = self.Example(self.U, self.dP, N, self.args)
        rhs = fd.make_rhs(ex)
        mat = fd.make_mat(ex)
        ps = np_solve(mat, rhs)
        pressure  = Pressure(ex, ps)
        flux = self.est_flux(ex, pressure.ps_2D)
        if plot:
            self.p_plot(ex, pressure, flux)
            
            velocity = Velocity(ex, ps)
            self.v_plot(ex, velocity, flux)
        return pressure, flux
    
    def fd_adj_solve(self, N, plot=True):
        ex = self.Example(self.U, self.dP, N, self.args)
        
        rhs = fd.make_rhs(ex)
        mat = fd.make_mat(ex)
        reyn_ps = np_solve(mat, rhs)
        adj_pressure = Adj_Pressure(ex, reyn_ps)
        adj_velocity = Adj_Velocity(ex, adj_pressure.ps_2D)
        flux = self.est_flux(ex, adj_pressure.ps_2D)
        
        if plot:
            self.p_plot(ex, adj_pressure, flux)
            self.v_plot(ex, adj_velocity, flux)
        
        #TODO: velocity for p_adj??
        
        return adj_pressure, flux

    def p_plot(self, ex, pressure, flux):
        paramstr = "$Re=0$, $Q=%.2f$, $U=%.1f$, $\Delta P=%.2f$"%(flux, self.U, self.dP)
        p_title = "" + paramstr
        p_labels = ["$p(x)$", "$x$","$y$"]
        graphics.plot_contour_mesh(pressure.ps_2D, ex.xs, ex.ys, p_title, p_labels, vmin=self.p_min, vmax=self.p_max, log_cmap=False)
    
    def v_plot(self, ex, velocity, flux):
        paramstr = "$Re=0$, $Q=%.2f$, $U=%.1f$, $\Delta P=%.2f$"%(flux, self.U, self.dP)
        v_title = "" + paramstr
        v_ax_labels =  ['$|(u,v)|_2$','$x$', '$y$'] 
        uv_mag = np.sqrt(velocity.vx**2 + velocity.vy**2)
        graphics.plot_stream_heat(velocity.vx, velocity.vy, ex.xs, ex.ys, uv_mag, v_title, v_ax_labels, vmin=0, vmax=self.vel_max, log_cmap=False)


def convg_rate(errs):
    n = len(errs)
    rates = np.zeros(n-1)
    
    for k in range(n-1):
        rates[k]=errs[k+1]/errs[k]
    
    return rates
    
def convg_pwl_fd(Example, U, dP, args, N0, dN, many,linthresh):
    inf_errs = np.zeros(many)
    l1_errs = np.zeros(many)
    l2_errs = np.zeros(many)
    Ns = np.zeros(many)
    N=N0
    solver = Reynolds_Solver(Example, U, dP, args)
    for k in range(many):
        Ns[k] = N
        pwl_p, flux = solver.pwl_solve(N, plot=False)
        fd_p, flux = solver.fd_solve(N, plot=False)
        size_k = len(pwl_p.ps_1D)
        
        p_err = np.abs(pwl_p.ps_1D - fd_p.ps_1D)
        inf_errs[k] = np.max(p_err)
        l1_errs[k] = np.sum(p_err)/size_k
        l2_errs[k] = np.sqrt(np.sum(p_err**2)/size_k)
        
        N*=dN
        
    title ='' #'Grid Convergence for Reynolds pressure \n Finite-Difference to Piecewise-Linear-Analytic'
    ax_labels=['$N$',  '$||p_{FD} - p_{PLA}||$']
    fun_labels=['$L_1$', '$L_2$', '$L_\infty$']

    graphics.plot_log_multi([l1_errs, l2_errs,inf_errs],Ns,title,fun_labels,ax_labels,linthresh, O1=2, O2=4)

    l1_rate = convg_rate(l1_errs)
    l2_rate = convg_rate(l2_errs)
    inf_rate = convg_rate(inf_errs)
    
    print("cnvg rates")
    print("l1: " + np.array2string(l1_rate))
    print("l2: " + np.array2string(l2_rate))
    print("linf" + np.array2string(inf_rate))

        
