# -*- coding: utf-8 -*-
"""
Created on Tue May 21 16:05:05 2024

@author: sarah
"""
import numpy as np
import graphics
import error
import reyn_pressures_pwlinear as pwl
from scipy.sparse.linalg import gmres as sp_gmres

import reyn_pressures_finDiff as fd
from numpy.linalg import solve as np_solve

from reyn_velocity import Velocity
from reyn_pressure import Pressure


class Reynolds_Solver: 
    def __init__(self, Example, U, dP):
        self.Example = Example
        
        self.U = U
        self.dP = dP
        
        self.vel_max = 4
        self.p_min=-70
        self.p_max=70
        
    def pwl_solve(self,N):
        ex = self.Example(self.U, self.dP, N)
        rhs = pwl.make_rhs(ex)
        linOp = pwl.pwlLinOp(ex)
        coefs, exit_code = sp_gmres(linOp, rhs, tol=1e-10)
        if exit_code != 0:
            print('gmres did not converge')
        ps, flux = pwl.make_ps(ex, coefs)
        pressure  = Pressure(ex, ps)
        velocity = Velocity(ex, ps)
        return pressure, flux, velocity
    
    def fd_solve(self, N):
        ex = self.Example(self.U, self.dP, N)
        rhs = fd.make_rhs(ex)
        mat = fd.make_mat(ex)
        ps = np_solve(mat, rhs)
        pressure  = Pressure(ex, ps)
        return pressure
    
    def solve_and_plot(self, N):
        pwl_pressure, flux, vel = self.pwl_solve(N)
        # fd_ps = self.fd_solve(N)
        ex = self.Example(self.U, self.dP, N)
        paramstr = "$Re=0$, $Q=%.2f$, $U=%.1f$, $\Delta P=%.2f$"%(flux, self.U, self.dP)
        p_title = "" + paramstr
        p_labels = ["$p(x)$", "$x$","$y$"]
        # graphics.plot_2D(pwl_pressure.ps_1D, ex.xs, p_title, p_labels, color='r')
        graphics.plot_contour_mesh(pwl_pressure.ps_2D, ex.xs, ex.ys, p_title, p_labels, vmin=self.p_min, vmax=self.p_max)
        # graphics.plot_2D(fd_ps, ex.xs, "fin-diff " + p_title, p_labels, color='r')
        
        v_title = "" + paramstr
        v_ax_labels =  ['$|(u,v)|_2$','$x$', '$y$'] 
        uv_mag = np.sqrt(vel.vx**2 + vel.vy**2)
        graphics.plot_stream_heat(vel.vx, vel.vy, ex.xs, ex.ys, uv_mag, v_title, v_ax_labels, vmin=0, vmax=self.vel_max)

def convg_pwl_fd(Example, U, dP, N0, dN, many):
    inf_errs = np.zeros(many)
    l1_errs = np.zeros(many)
    l2_errs = np.zeros(many)
    Ns = np.zeros(many)
    N=N0
    
    solver = Reynolds_Solver(Example, U, dP)
    for k in range(many):
        Ns[k] = N
        pwl_p, flux, vel = solver.pwl_solve(N)
        fd_p = solver.fd_solve(N)
        
        errs = np.abs(pwl_p.ps_1D - fd_p.ps_1D)
        inf_errs[k] = np.max(errs)
        l1_errs[k] = np.sum(errs)
        l2_errs[k] = np.sqrt(np.sum(errs**2))
        N*=2
        
    title = 'Grid Convergence for Reynolds pressure \n Finite-Difference to Piecewise-Linear-Analytic'
    ax_labels=['$N$',  '$||p_{FD} - p_{PWL}||$']
    fun_labels=['$L_1$', '$L_2$', '$L_\infty$']
    o1=2
    o2=4
    graphics.plot_log_multi([l1_errs, l2_errs,inf_errs],Ns,title,fun_labels,ax_labels,linthresh=1e-5, O1=o1, O2=o2)
    
    l1_rate = error.convg_rate(l1_errs)
    l2_rate = error.convg_rate(l2_errs)
    inf_rate = error.convg_rate(inf_errs)
    
    print("cnvg rates")
    print("l1: " + np.array2string(l1_rate))
    print("l2: " + np.array2string(l2_rate))
    print("linf" + np.array2string(inf_rate))

        