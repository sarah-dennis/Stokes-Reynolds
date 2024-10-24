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

from reyn_velocity import Velocity

class Reynolds_Solver: 
    def __init__(self, example, U, dP):
        self.example = example
        
        self.U = U
        self.dP = dP
        
    def pwl_solve(self,N):
        ex = self.example(self.U, self.dP, N)
        rhs = pwl.make_rhs(ex)
        linOp = pwl.pwlLinOp(ex)
        coefs, exit_code = sp_gmres(linOp, rhs, tol=1e-16)
        ps, flux = pwl.make_ps(ex, coefs)
        vel = Velocity(ex, ps)
        return ps, flux, vel
    
    def fd_solve(self, N):
        ex = self.example(self.U, self.dP, N)
        rhs = fd.make_rhs(ex)
        mat = fd.make_mat(ex)
        ps = np_solve(mat, rhs)
        return ps
    
    def solve_plot(self, N):
        pwl_ps, flux, vel = self.pwl_solve(N)
        # fd_ps = self.fd_solve(N)
        ex = self.example(self.U, self.dP, N)
        paramstr = "$Re=0$, $Q=%.2f$, $U=%.1f$, $\Delta P=%.2f$"%(flux, self.U, self.dP)
        p_title = "Pressure $p(x)$: \n" + paramstr
        p_labels = ["$x$", "Pressure $p(x)$"]
        graphics.plot_2D(pwl_ps, ex.xs, p_title, p_labels, color='r')
        # graphics.plot_2D(fd_ps, ex.xs, "Fin-Diff " + p_title, p_labels, color='r')
        
        #y-axis reverse
        vy=np.flip(vel.vy, 0)
        vx=np.flip(vel.vx, 0)
        v_title = "Velocity $(u,v)$: \n" + paramstr
        v_ax_labels =  ['$|(u,v)|_2$','$x$', '$y$'] 
        uv_mag = np.sqrt(vx**2 + vy**2)
        vmax = 2.5*np.median(uv_mag) #TODO explodes at vertical edges
        graphics.plot_stream_heat(vx, vy, ex.xs, ex.ys, uv_mag, v_title, v_ax_labels, vmin=0, vmax=vmax)

    def cnvg(self, N0, dN, many):
        
        inf_errs = np.zeros(many)
        l1_errs = np.zeros(many)
        l2_errs = np.zeros(many)
        Ns = np.zeros(many)
        N=N0
        for k in range(many):
            Ns[k] = N
            pwl_ps, flux, vel = self.pwl_solve(N)
            fd_ps = self.fd_solve(N)
            errs = np.abs(pwl_ps - fd_ps)
            inf_errs[k] = np.max(errs)
            l1_errs[k] = np.sum(errs)
            l2_errs[k] = np.sqrt(np.sum(errs**2))
            N*=2
        title = 'error FD to PWL Reynolds'
        ax_labels=['N',  'err']
        fun_labels=['l1', 'l2', 'inf']
        graphics.plot_log_multi([l1_errs, l2_errs,inf_errs],Ns,title,fun_labels,ax_labels,linthresh=1e-8)
        
        l1_rate = convg_rate(l1_errs)
        l2_rate = convg_rate(l2_errs)
        inf_rate = convg_rate(inf_errs)
        
        print("cnvg rates")
        print("l1: " + np.array2string(l1_rate))
        print("l2: " + np.array2string(l2_rate))
        print("linf" + np.array2string(inf_rate))


            
def convg_rate(errs):
    n = len(errs)
    rates = np.zeros(n-1)
    for k in range(n-1):
        rates[k]=errs[k+1]/errs[k]
    
    return rates

        