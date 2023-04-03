# -*- coding: utf-8 -*-
"""
Created on Tue Mar  7 16:16:13 2023

@author: sarah
"""
import numpy as np
import Reynolds_1D as ry
import domain as dfd
import examples_1D as eg
from matplotlib import pyplot as pp
import _graphics as g
import exactPressures as ps
import heights as hg
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
    
    pp.figure(trials+1)

    pp.loglog(dxs, dxs, color='b', label='$\mathcal{O}(dx)$')
    pp.loglog(dxs, dxs_sqr, color='g', label='$\mathcal{O}(dx^2)$')
    pp.loglog(dxs, infNorms, color='r', label='$L_{\infty}$ Norm Error')
    
    pp.xlabel('Grid spacing $dx$')

    pp.legend(loc='lower right')
    pp.title("Convergence for %s"%(height_k.h_str), fontweight ="bold")


#------------------------------------------------------------------------------
# Parameter variation
#------------------------------------------------------------------------------
#square wave: period --> 0
def vary_n_steps(trials=20, n_steps_0=5):
    # ry.domain <- (x0, xf, Nx, BC, U, eta, dx)
    n_steps_k = n_steps_0
    r = 0.1

    v = np.zeros(trials)
    x = np.zeros(trials)
    
    for k in range(trials):
        height_k, pressure_k = eg.squareWave(ry.domain, ry.p0, ry.pN, n_steps_k, r)
        
        
        
        #err_k, ps_k = ry.solve(ry.domain, height_k, pressure_k, 0, 1)
        
        #v[k] = np.abs(ps_k[1]-ps_k[0])
        #title = "Pressure Drop as steps width decreases"
        #y_axis = "$|p_1-p_0|$"
        
        #v[k] = np.max(ps_k)
        
        
        v[k] = np.max(pressure_k.ps)-ry.p0
        
        
        title = "Max Pressure as step width decreases"
        y_axis = "$p_{max}$"
        
        x[k] = n_steps_k
                
        n_steps_k = n_steps_k * 2 + 1
        
    #g.plot_2D(v, x, title, y_axis, "number of steps over $x\in[%.1f, %.1f]$"%(ry.x0, ry.xf))
    pp.loglog( v, x, color='b', label='$\mathcal{O}(dx)$')


def vary_n_steps_condNum(trials=4, n_steps_0=5):
    n_steps_k = n_steps_0
    r = 0.1
    h_avg = 0.2
         
    v = np.zeros(trials)
    x = np.zeros(trials)
    
    for k in range(trials):
        sqr_height_k = hg.SquareWaveHeight(ry.domain, h_avg, r, n_steps_k)
        sqr_pressure_k = ps.SquareWavePressure(ry.domain, sqr_height_k, ry.p0, ry.pN)
        v[k] = np.linalg.cond(sqr_pressure_k.M)
        x[k] = n_steps_k
        
        title = "Condition Number as step width decreases"
        y_axis = "Condition Number"
        
        x[k] = n_steps_k
                
        n_steps_k = n_steps_k * 2 + 1
        
    g.plot_2D(v, x, title, y_axis, "number of steps over $x\in[%.1f, %.1f]$"%(ry.x0, ry.xf))


#square wave: radous --> 0
def vary_r(trials=15, r_0=0.05):
    # ry.domain <- (x0, xf, Nx, BC, U, eta, dx)
    n_steps = 5
    r_k = r_0

    v = np.zeros(trials)
    x = np.zeros(trials)
    
    for k in range(trials):
        height_k, pressure_k = eg.squareWave(ry.domain, ry.p0, ry.pN, n_steps, r_k)
        err_k, ps_k = ry.solve(ry.domain, height_k, pressure_k, 0, 1)
        
        #v[k] = np.abs(ps_k[1]-ps_k[0])
        #title = "Pressure Drop as steps width decreases"
        #y_axis = "$|p_1-p_0|$"
        
        v[k] = np.max(ps_k)
        title = "Max Pressure as step radius decreases"
        y_axis = "$p_{max}$"
        
        x[k] = r_k
                
        r_k /= 2
        
    g.plot_2D(v, x, title, y_axis, "height radius")






