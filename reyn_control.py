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

class Reynolds_Solver: 
    def __init__(self, Example, U, dP):
        self.Example = Example
        
        self.U = U
        self.dP = dP
        
    def pwl_solve(self,N):
        ex = self.Example(self.U, self.dP, N)
        rhs = pwl.make_rhs(ex)
        linOp = pwl.pwlLinOp(ex)
        coefs, exit_code = sp_gmres(linOp, rhs, tol=1e-10)
        if exit_code != 0:
            print('gmres did not converge')
        ps, flux = pwl.make_ps(ex, coefs)
        vel = Velocity(ex, ps)
        return ps, flux, vel
    
    def fd_solve(self, N):
        ex = self.Example(self.U, self.dP, N)
        rhs = fd.make_rhs(ex)
        mat = fd.make_mat(ex)
        ps = np_solve(mat, rhs)
        return ps
    
    def solve_and_plot(self, N):
        pwl_ps, flux, vel = self.pwl_solve(N)
        # fd_ps = self.fd_solve(N)
        ex = self.Example(self.U, self.dP, N)
        paramstr = "$Re=0$, $Q=%.2f$, $U=%.1f$, $\Delta P=%.2f$"%(flux, self.U, self.dP)
        p_title = "Pressure $p(x)$: \n" + paramstr
        p_labels = ["$x$", "Pressure $p(x)$"]
        graphics.plot_2D(pwl_ps, ex.xs, p_title, p_labels, color='r')
        # graphics.plot_2D(fd_ps, ex.xs, "fin-diff " + p_title, p_labels, color='r')
        
        #y-axis reverse
        vy=np.flip(vel.vy, 0)
        vx=np.flip(vel.vx, 0)
        v_title = "Velocity $(u,v)$: \n" + paramstr
        v_ax_labels =  ['$|(u,v)|_2$','$x$', '$y$'] 
        uv_mag = np.sqrt(vx**2 + vy**2)
        vmax = 2.5*np.median(uv_mag) #TODO velocity explodes at vertical edges
        graphics.plot_stream_heat(vx, vy, ex.xs, ex.ys, uv_mag, v_title, v_ax_labels, vmin=0, vmax=vmax)

def convg_pwl_fd(Example, U, dP, N0, dN, many):
    inf_errs = np.zeros(many)
    l1_errs = np.zeros(many)
    l2_errs = np.zeros(many)
    Ns = np.zeros(many)
    N=N0
    
    solver = Reynolds_Solver(Example, U, dP)
    for k in range(many):
        Ns[k] = N
        pwl_ps, flux, vel = solver.pwl_solve(N)
        fd_ps = solver.fd_solve(N)
        
        errs = np.abs(pwl_ps - fd_ps)
        inf_errs[k] = np.max(errs)
        l1_errs[k] = np.sum(errs)
        l2_errs[k] = np.sqrt(np.sum(errs**2))
        N*=2
    title = 'error FD to PWL Reynolds'
    ax_labels=['N',  'err']
    fun_labels=['l1', 'l2', 'inf']
    graphics.plot_log_multi([l1_errs, l2_errs,inf_errs],Ns,title,fun_labels,ax_labels,linthresh=1e-8)
    
    l1_rate = error.convg_rate(l1_errs)
    l2_rate = error.convg_rate(l2_errs)
    inf_rate = error.convg_rate(inf_errs)
    
    print("cnvg rates")
    print("l1: " + np.array2string(l1_rate))
    print("l2: " + np.array2string(l2_rate))
    print("linf" + np.array2string(inf_rate))


        