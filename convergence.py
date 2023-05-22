# -*- coding: utf-8 -*-
"""
Created on Tue Mar  7 16:16:13 2023

@author: sarah
"""
import numpy as np
import Reynolds_1D as ry
import domain as dfd
import examples_1D as eg
import _graphics as g

#------------------------------------------------------------------------------
# Convergence
#------------------------------------------------------------------------------   
#BCs = [0: periodic, 1: fixed]
#RHS = [0: reynolds, 1: manufactured]

# dx -> 0
def conveg(trials=10, N0=50, BC="fixed", RHS=0):
    Nx_k = N0
    infNorms = np.zeros(trials)
    dxs = np.zeros(trials)
    dxs_sqr = np.zeros(trials)
    fig=0
    for k in range(trials):
        if k == trials-1: fig=1
        
        domain_k = dfd.Domain(ry.x0, ry.xf, ry.eta, ry.U, Nx_k, ry.BC)
        
        #---- EXAMPLES -----------------------------------------
        #height_k, pressure_k = eg.wedge(domain_k)
        #height_k, pressure_k = eg.corrugated(domain_k)
        #height_k, pressure_k = eg.step(domain_k)
        #height_k, pressure_k = eg.twoStep(domain_k)
        height_k, pressure_k = eg.squareWave(domain_k, ry.p0, ry.pN)
        #-------------------------------------------------------
        err_k, ps_k = ry.solve(domain_k, height_k, pressure_k, RHS, FIG=fig)
        
        infNorms[k] = err_k
        
        dxs[k] = domain_k.dx
        dxs_sqr[k] = dxs[k]**2
        
        Nx_k *= 2
    
    
    labels = ['$\mathcal{O}(dx)$', '$\mathcal{O}(dx^2)$', '$L_{\infty}$ Norm Error']
    fs = [dxs, dxs_sqr, infNorms]
    x_axis = '$dx$'
    title = "Convergence for %s"%(height_k.h_str)
    g.plot_log_multi(fs, dxs, title, labels, x_axis)


#------------------------------------------------------------------------------
# Parameter variation for square wave
#------------------------------------------------------------------------------

def vary_nSteps_pMax(trials=7, n_steps_0=5):
    # ry.domain <- (x0, xf, Nx, BC, U, eta, dx)
    n_steps_k = n_steps_0
    r = 0.1

    v = np.zeros(trials)
    x = np.zeros(trials)
    
    for k in range(trials):
        height_k, pressure_k = eg.squareWave(ry.domain, ry.p0, ry.pN, n_steps_k, r)
        
        #err_k, pressure_num_k = ry.solve(ry.domain, height_k, pressure_k, 0, 1)
        
        v[k] = np.max(pressure_k.ps)-ry.p0

        x[k] = n_steps_k
                
        n_steps_k = n_steps_k * 2 + 1
    title = "Max Pressure as step width decreases"
    y_axis = "$p_{max}$"
    x_axis = "$dx$"  
    g.plot_log(v, x, title,  x_axis, y_axis)


def vary_nSteps_condNum(trials=6, n_steps_0=5):
    n_steps_k = n_steps_0
    r = 0.1
    h_avg = 0.2
         
    v = np.zeros(trials)
    u = np.zeros(trials)
    x = np.zeros(trials)
    
    for k in range(trials):
        
        height_k, pressure_k = eg.squareWave(ry.domain, ry.p0, ry.pN, n_steps_k, r, h_avg)
        
        u[k] = np.linalg.cond(pressure_k.M)
        v[k] = np.linalg.cond(pressure_k.M_inv)
        x[k] = n_steps_k
    
        n_steps_k = n_steps_k * 2 + 1
        
    title = "Condition Number vs Matrix Size"
    y_axis = "Condition Number"
    x_axis = "n steps for $x\in[%.1f, %.1f]$"%(ry.x0, ry.xf)
    labels = ["M", "M_inv"]
    g.plot_2D_multi([u,v], x, title, labels, [x_axis, y_axis] )


def vary_r_pMax(trials=10, r_0=0.05):
    # ry.domain <- (x0, xf, Nx, BC, U, eta, dx)
    n_steps = 5
    r_k = r_0

    v = np.zeros(trials)
    x = np.zeros(trials)
    
    for k in range(trials):
        height_k, pressure_k = eg.squareWave(ry.domain, ry.p0, ry.pN, n_steps, r_k)
        
        #err_k, ps_k = ry.solve(ry.domain, height_k, pressure_k, 0, 1)
        
        v[k] = np.max(pressure_k.ps) - ry.p0

        x[k] = r_k
                
        r_k /= 2
        
    title = "Max Pressure as step radius decreases"
    y_axis = "$p_{max}$"
    g.plot_2D(v, x, title, y_axis, "height radius")






