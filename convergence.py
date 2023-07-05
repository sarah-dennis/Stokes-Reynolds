# -*- coding: utf-8 -*-
"""
Created on Tue Mar  7 16:16:13 2023

@author: sarah
"""
import numpy as np

import domain as dfd
import Reynolds_1D as ry
import examples_1D as eg
import _graphics as g
import groundControl as gc


#------------------------------------------------------------------------------
# Convergence of Numerical Reynolds to Analytic
#------------------------------------------------------------------------------   
#BCs = [0: periodic, 1: fixed]
#RHS = [0: reynolds, 1: manufactured]

# Nx -> large, dx -> 0
def vary_Nx_numErr(trials=6, N0 = 50):
    Nx_k = N0
    infNorms = np.zeros(trials)
    dxs = np.zeros(trials)
    dxs_sqr = np.zeros(trials)

    for k in range(trials):
        
        domain_k = dfd.Domain(gc.x0, gc.xf, gc.eta, gc.U, Nx_k, gc.BC)
        
        #---- EXAMPLES -----------------------------------------
        #height_k, pressure_k = eg.wedge(domain_k)
        #height_k, pressure_k = eg.corrugated(domain_k)
        #height_k, pressure_k = eg.step(domain_k)
        #height_k, pressure_k = eg.twoStep(domain_k)
        height_k, pressure_k = eg.squareWave(domain_k, gc.p0, gc.pN)
        #-------------------------------------------------------
        numPressure_k = ry.solve(domain_k, height_k, gc.p0, gc.pN)
        
        infNorms[k] = np.max(np.abs(pressure_k.ps - numPressure_k))
        
        dxs[k] = domain_k.dx
        dxs_sqr[k] = dxs[k]**2
        
        Nx_k *= 2
    
    labels = ['$\mathcal{O}(dx)$', '$\mathcal{O}(dx^2)$', '$L_{\infty}$ Norm Error']
    fs = [dxs, dxs_sqr, infNorms]
    x_axis = '$dx$'
    y_axis = 'Error'
    title = "Convergence for %s"%(height_k.h_str)
    g.plot_log_multi(fs, dxs, title, labels, [x_axis, y_axis])


#------------------------------------------------------------------------------
# Parameter variation for square wave
#------------------------------------------------------------------------------

def vary_nSteps_pMax(trials=7, n_steps_0=5):
    # ry.domain <- (x0, xf, Nx, BC, U, eta, dx)
    n_steps_k = n_steps_0
    r = 0.1
    h_avg = 0.2

    v = np.zeros(trials)
    x = np.zeros(trials)
    
    for k in range(trials):
        height_k, pressure_k = eg.squareWave(gc.domain, gc.p0, gc.pN, n_steps_k, r, h_avg)
        
        #err_k, pressure_num_k = ry.solve(ry.domain, height_k, pressure_k, 0, 1)
        
        v[k] = np.max(pressure_k.ps)-gc.p0

        x[k] = n_steps_k
                
        n_steps_k = n_steps_k * 2 + 1
    title = "Max Pressure as step width decreases"
    y_axis = "$p_{max}$"
    x_axis = "$dx$"  
    g.plot_log(v, x, title,  x_axis, y_axis)
    
    

def vary_nSteps_time(trials=9, n_steps_0=101):
    n_steps_k = n_steps_0
    r = 0.001
    h_avg = 0.1

    u = np.zeros(trials)
    w = np.zeros(trials)
    x = np.zeros(trials)
    v = np.zeros(trials)
    
    for k in range(trials):
        
        height_k, pressure_py_k = eg.squareWave_pySolve(gc.domain, gc.p0, gc.pN, n_steps_k, r, h_avg)
        height_k, pressure_lu_k = eg.squareWave_schurLUSolve(gc.domain, gc.p0, gc.pN, n_steps_k, r, h_avg)
        
        u[k] = pressure_lu_k.time
        w[k] = pressure_py_k.time
        v[k] = n_steps_k**2
        x[k] = n_steps_k
    
        n_steps_k = int(n_steps_k * 2+1)
        
    title = "Solve Time vs Matrix Size"
    y_axis = "Time"
    x_axis = "n steps for $x\in[%.1f, %.1f]$"%(gc.x0, gc.xf)
    g.plot_log_multi([u,w, v], x, title, ["LU runtime", "python runtime", "$\mathcal{O}(n^2)$"], [x_axis, y_axis])
    


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






