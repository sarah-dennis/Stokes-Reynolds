# -*- coding: utf-8 -*-
"""
Created on Tue Mar  7 16:16:13 2023

@author: sarah
"""
import numpy as np

import domain as dfd
import Reynolds_1D as ry
import examples_1D as eg
import heights as hgt
import pressures as prs
import _graphics as g
import control as gc
import csv


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
    
    

def vary_nSteps_time(trials=1, repeats=1, n_steps_0=65001):
    bigO_const = 10**-6
    n_steps_k = n_steps_0
    r = 0.001
    h_avg = 0.1
    
    for rep in range(repeats):
        print("\n Round %d of %d"%(rep+1, repeats))
        n_steps_k = n_steps_0
        
        lu = np.zeros(trials)
        
        # inv = np.zeros(trials)
        py = np.zeros(trials)
        
        n = np.zeros(trials)
        nsqr = np.zeros(trials)

        py_fail = False
        lu_fail = False
        
        for k in range(trials):
            print("\n Trial %d of %d: N steps = %d"%(k+1, trials, n_steps_k))
            height_k = hgt.SquareWaveHeight(gc.domain, h_avg, r, n_steps_k)
            
            # if not py_fail:
            #     try:
            #         pressure_py_k = prs.SquareWavePressure_pySolve(gc.domain, height_k, gc.p0, gc.pN)
            #         py[k] = pressure_py_k.time
            #     except MemoryError as error:
            #         print("In Python solve: ", error)
            #         py_fail = True
                    
            if not lu_fail:
                try:
                    pressure_lu_k = prs.SquareWavePressure_schurLUSolve(gc.domain, height_k, gc.p0, gc.pN)
                    lu[k] = pressure_lu_k.time
                except MemoryError as error:
                    print("In LU solve: ", error)
                    lu_fail = True
            
                
            # pressure_inv_k = prs.SquareWavePressure_schurInvSolve(gc.domain, height_k, gc.p0, gc.pN)
            # inv[k] = pressure_inv_k.time
            
            n[k] = n_steps_k
            nsqr[k] = bigO_const * n_steps_k**2
            n_steps_k = int(n_steps_k * 2+1)
    
        file_name = "sw_solveTimes_%s.csv"%(rep+1)
        np.savetxt(file_name, np.array([n, lu]).T, delimiter =", ", fmt='%.5f')
        
        plot_title = "Solve Time vs Matrix Size (Trial %d)"%(rep + 1)
        y_axis = "Time"
        x_axis = "$N_{steps}$"
        # funs = [nsqr, lu, py]
        funs = [nsqr, lu]
        # fun_labels = ["$\mathcal{O}(n^2)$","LU solve time", "numpy solve time"]
        fun_labels = ["$\mathcal{O}(n^2)$","LU solve time"]
        # g.plot_log_multi([nsqr, lu, py, inv], n, plot_title, ["$\mathcal{O}(n^2)$","LU runtime", "python runtime", "SchurInv runtime",], [x_axis, y_axis])
        g.plot_log_multi(funs, n, plot_title, fun_labels, [x_axis, y_axis])

        


def vary_r_pMax(trials=10, r_0=0.05):
    # ry.domain <- (x0, xf, Nx, BC, U, eta, dx)
    n_steps = 5
    r_k = r_0

    v = np.zeros(trials)
    x = np.zeros(trials)
    
    for k in range(trials):
        height_k, pressure_k = eg.squareWave(gc.domain, gc.p0, gc.pN, n_steps, r_k)
        
        #err_k, ps_k = ry.solve(ry.domain, height_k, pressure_k, 0, 1)
        
        v[k] = np.max(pressure_k.ps) - ry.p0

        x[k] = r_k
                
        r_k /= 2
        
    title = "Max Pressure as step radius decreases"
    y_axis = "$p_{max}$"
    g.plot_2D(v, x, title, y_axis, "height radius")






