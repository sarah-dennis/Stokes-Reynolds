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


Nx = 5000   # Number of Grid points

BC = "fixed" # Boundary Condition in x (alt. "periodic")

domain = dm.Domain(x0, xf, eta, U, Nx, BC)


#------------------------------------------------------------------------------
# Height & Pressure
#------------------------------------------------------------------------------
p0 = 0
pN = 0

n_steps = 1527
# Overflow when finding (thetas, phis) around 1560 steps

#---------------------------------------------------------------------------
# Initialize Height and Pressure

#height, pressure = flat(domain, p0, pN)
#height, pressure = wedge(domain, p0, pN)
#height, pressure = corrugated(domain, p0, pN)
#height, pressure = step(domain, p0, pN)
#height, pressure = twoStep(domain, p0, pN)

height, pressure_inv = egs.squareWave_schurInvSolve(domain, p0, pN, n_steps)
height, pressure_py = egs.squareWave_pySolve(domain, p0, pN, n_steps)
height, pressure_lu = egs.squareWave_schurLUSolve(domain, p0, pN, n_steps)

#---------------------------------------------------------------------------
# Plotting 

p_h_title = "Analytic Pressure for %s"%height.h_str
p_h_labels = ["Pressure $p(x)$", "Height $h(x)$", "$x$"]
graph.plot_2D_twin(pressure_py.ps, height.hs, domain.xs, p_h_title, p_h_labels)

p_p_title = "LU Solve vs Python Solve of Pressure for %s"%height.h_str
p_p_labels = ["LU Pressure", "Python Pressure", "x"]
graph.plot_2D_twin(pressure_lu.ps, pressure_py.ps, domain.xs, p_p_title, p_p_labels)

p_p_title = "SchurInv Solve vs PySolve of Pressure for %s"%height.h_str
p_p_labels = ["LU Pressure", "Python Pressure", "x"]
graph.plot_2D_twin(pressure_inv.ps, pressure_py.ps, domain.xs, p_p_title, p_p_labels)

#---------------------------------------------------------------------------
# Numerical Solution

# num_pressure = ry.solve(domain, height, p0, pN)

# num_err = np.max(np.abs(num_pressure - pressure.ps))
# print("Analytic to Numerical Error: %.3f"%num_err)

# nump_h_title = "Reynolds Numerical Pressure for %s"%height.h_str
# graph.plot_2D_twin(num_pressure.ps, height.hs, domain.xs, nump_h_title)




