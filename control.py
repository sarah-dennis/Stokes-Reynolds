#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 22 15:05:40 2023

@author: sarahdennis
"""

import numpy as np

import domain as dm
import Reynolds_1D as ry
import examples_1D as ex
import _graphics as graph


#------------------------------------------------------------------------------
# Domain & Discretization
#------------------------------------------------------------------------------
x0 = 0      # left boundary 
xf = 1      # right boundary
U = 2     # surface velocity
eta = 1     # viscosity


Nx = 1000   # Number of Grid points

BC = "fixed" # Boundary Condition in x (alt. "periodic")

domain = dm.Domain(x0, xf, eta, U, Nx, BC)

n_steps = 3

#------------------------------------------------------------------------------
# Height & Pressure
#-------------------------------------------------------------------------------
p0 = 0
pN = 0

#---------------------------------------------------------------------------
# Analytic Solution

# height, pressure = ex.flat(domain, p0, pN)
# height, pressure = ex.wedge(domain, p0, pN)
height, pressure = ex.corrugated(domain, p0, pN)
# height, pressure = ex.step(domain, p0, pN)
# height, pressure = ex.twoStep(domain, p0, pN)

# height, pressure = ex.squareWave_schurInvSolve(domain, p0, pN, n_steps)
# height, pressure = ex.squareWave_schurLUSolve(domain, p0, pN, n_steps)


ys = np.linspace(0, 1.2*height.h_max, domain.Nx)

v_x = 1/(2*domain.eta) * pressure.pxs * (ys**2 - ys*height.h_max) + ys*domain.U/height.h_max
v_y = domain.U * height.hxs

phv_title = "Pressure and Velocity for %s"%height.h_str
phv_fun_labels = ['velocity $(v_x, v_y)$', 'pressure $p(x)$', 'height $h(x)$']
phv_ax_labels =  ['x', 'y']

graph.plot_phv(pressure.ps, height.hs, v_x, v_y, domain.xs, ys, phv_title, phv_fun_labels,  phv_ax_labels)

    
#---------------------------------------------------------------------------
# Numerical Solution

# height, pressure_gmres = ex.squareWave_gmresSolve(domain, p0, pN, n_steps)
# height, pressure_py = ex.squareWave_pySolve(domain, p0, pN, n_steps)

# pressure_ry = ry.solve(domain, height, p0, pN)

# num_err = np.max(np.abs(num_pressure - pressure.ps))
# print("Analytic to Numerical Error: %.3f"%num_err)

# nump_h_title = "Reynolds Numerical Pressure for %s"%height.h_str
# graph.plot_2D_twin(num_pressure.ps, height.hs, domain.xs, nump_h_title)


#------------------------------------------------------------------------------
# Error

# err = la.norm(pressure_lu.ps - pressure_gmres.ps)
# print("analytic to numerical error %.5f"%(err))

#-------------------------
# fluid velocity



#---------------------------------------------------------------------------
# Plotting 


p_h_title = "Pressure (%s) and Height for %s"%(pressure.p_str, height.h_str)
p_h_labels = ["Pressure $p(x)$", "Height $h(x)$", "$x$"]
graph.plot_2D_twin(pressure.ps, height.hs, domain.xs, p_h_title, p_h_labels)

# graph.plot_2D_twin(pressure_lu.ps, height.hs, domain.xs, p_h_title, p_h_labels)

# graph.plot_2D_twin(pressure_lu.ps - pressure_gmres.ps, height.hs, domain.xs, "LU - Gmres Pressure", p_h_labels)




