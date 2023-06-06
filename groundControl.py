#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 22 15:05:40 2023

@author: sarahdennis
"""

import domain as dm
import pressures as prs
import heights as hgt
import numpy as np
import _graphics as graph

#------------------------------------------------------------------------------
# Domain & Discretization
#------------------------------------------------------------------------------
x0 = 0      # left boundary 
xf = 1      # right boundary
U = 1       # surface velocity
eta = 1     # viscosity


Nx = 10000   # Number of Grid points

BC = "fixed" # Boundary Condition in x (alt. "periodic")

domain = dm.Domain(x0, xf, eta, U, Nx, BC)


#------------------------------------------------------------------------------
# Height & Pressure
#------------------------------------------------------------------------------
p0 = 0
pN = 0

n_steps = 2001


# Initialise Height and Pressure

#height, pressure = flat(domain, p0, pN)
#height, pressure = wedge(domain, p0, pN)
#height, pressure = corrugated(domain, p0, pN)
#height, pressure = step(domain, p0, pN)
#height, pressure = twoStep(domain, p0, pN)

# ------ Paste from examples_1D here ------------------------------------------
def squareWave(domain, p0, pN, n_steps=7 , r=0.05, h_avg=0.1):
    if domain.Nx < n_steps * 3:
        print("Warning: Nx < nsteps * 3")
        
    print("Loading %d-step Square Wave \n"%(n_steps))
    height = hgt.SquareWaveHeight(domain, h_avg, r, n_steps)
    pressure = prs.SquareWavePressure(domain, height, p0, pN)

    return height, pressure

height, pressure = squareWave(domain, p0, pN, n_steps)


p_h_title = "Reynolds Analytic Pressure for %s"%height.h_str
graph.plot_2D_twin(pressure.ps, height.hs, domain.xs, p_h_title)
# ------------------------------------------------------------------------------

num_pressure = prs.ReynoldsPressure(domain, height, p0, pN)
nump_h_title = "Reynolds Numerical Pressure for %s"%height.h_str
graph.plot_2D_twin(num_pressure.ps, height.hs, domain.xs, nump_h_title)

max_err = np.max(np.abs(num_pressure.ps - pressure.ps))

# print("Analytic to Numerical Error: %.3f"%max_err)
# height.plot(domain)
