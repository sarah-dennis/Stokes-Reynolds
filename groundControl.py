#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 22 15:05:40 2023

@author: sarahdennis
"""

import domain as dm
import Reynolds_1D as ry
import numpy as np
import examples_1D as egs
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

n_steps = 3001
# Overflow when finding (thetas, phis) around 1560 steps

#---------------------------------------------------------------------------
# Initialize Height and Pressure

#height, pressure = flat(domain, p0, pN)
#height, pressure = wedge(domain, p0, pN)
#height, pressure = corrugated(domain, p0, pN)
#height, pressure = step(domain, p0, pN)
#height, pressure = twoStep(domain, p0, pN)

# height, pressure = egs.squareWave_schurInvSolv(domain, p0, pN, n_steps)
height, pressure_py = egs.squareWave_pySolve(domain, p0, pN, n_steps)

# anyl_err = np.max(np.abs(pressure.ps - pressure_py.ps))
# print("SchurComp Solve to Python Solve Error: %.3f \n"%anyl_err)

#---------------------------------------------------------------------------
# Plotting 
# height.plot(domain)
# pressure.plot(domain)

# p_h_title = "Reynolds Analytic Pressure for %s"%height.h_str
# graph.plot_2D_twin(pressure.ps, height.hs, domain.xs, p_h_title)

#---------------------------------------------------------------------------
# Numerical Solution

# num_pressure = ry.solve(domain, height, p0, pN)

# num_err = np.max(np.abs(num_pressure - pressure.ps))
# print("Analytic to Numerical Error: %.3f"%num_err)

# nump_h_title = "Reynolds Numerical Pressure for %s"%height.h_str
# graph.plot_2D_twin(num_pressure.ps, height.hs, domain.xs, nump_h_title)




