#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 22 15:05:40 2023

@author: sarahdennis
"""

import domain as dm
import examples_1D as egs
import pressures as prs
import numpy as np
import _graphics as graph

#------------------------------------------------------------------------------
# Domain & Discretization
#------------------------------------------------------------------------------
x0 = 0      # left boundary 
xf = 1      # right boundary
U = 1       # surface velocity
eta = 1     # viscosity


Nx = 4000   # Number of Grid points

BC = "fixed" # Boundary Condition in x (alt. "periodic")

domain = dm.Domain(x0, xf, eta, U, Nx, BC)


#------------------------------------------------------------------------------
# Height & Pressure
#------------------------------------------------------------------------------
p0 = 0
pN = 0


# Initialise Height and Pressure

#height, pressure = egs.flat(domain, p0, pN)
#height, pressure = egs.wedge(domain, p0, pN)
#height, pressure = egs.corrugated(domain, p0, pN)
#height, pressure = egs.step(domain, p0, pN)
#height, pressure = egs.twoStep(domain, p0, pN)
height, pressure = egs.squareWave(domain, p0, pN)

p_h_title = "Reynolds Analytic Pressure for %s"%height.h_str
graph.plot_2D_twin(pressure.ps, height.hs, domain.xs, p_h_title)

# num_pressure = prs.ReynoldsPressure(domain, height, p0, pN)
# nump_h_title = "Reynolds Numerical Pressure for %s"%height.h_str
# graph.plot_2D_twin(num_pressure.ps, height.hs, domain.xs, nump_h_title)

# max_err = np.max(np.abs(num_pressure.ps - pressure.ps))

# print("Analytic to Numerical Error: %.3f"%max_err)

