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
#------------------------------------------------------------------------------
# Convergence
#------------------------------------------------------------------------------   
#BCs = [0: periodic, 1: fixed]
#RHS = [0: reynolds, 1: manufactured]

def conveg(trials=10, N0=10, BC=1, RHS=0):
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
        height_k, pressure_k = eg.squareWave(domain_k)
        #-------------------------------------------------------

        
        infNorms[k] = ry.solve(domain_k, height_k, pressure_k, RHS, FIG=fig)
        
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
    