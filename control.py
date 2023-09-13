#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 22 15:05:40 2023

@author: sarahdennis
"""

import numpy as np

import cProfile as cpr

import domain as dm
import Reynolds_1D as ry
import examples_1D as ex
import _graphics as graph


#------------------------------------------------------------------------------
# Domain & Discretization
#------------------------------------------------------------------------------
x0 = 0      # left boundary 
xf = 1      # right boundary
U = 5     # surface velocity
eta = 1     # viscosity


Nx = 500   # Number of Grid points

BC = "fixed" # Boundary Condition in x (alt. "periodic")

domain = dm.Domain(x0, xf, eta, U, Nx, BC)

n_steps = 11

#------------------------------------------------------------------------------
# Height & Pressure
#-------------------------------------------------------------------------------
p0 = 0
pN = 0

#---------------------------------------------------------------------------
# Analytic Solution

# height, pressure = ex.flat(domain, p0, pN)
# height, pressure = ex.wedge(domain, p0, pN)
# height, pressure = ex.corrugated(domain, p0, pN)
# height, pressure = ex.step(domain, p0, pN)
# height, pressure = ex.twoStep(domain, p0, pN)

# height, pressure = ex.squareWave_schurInvSolve(domain, p0, pN, n_steps)
height, pressure = ex.squareWave_schurLUSolve(domain, p0, pN, n_steps)


#---------------------------------------------------------------------------
# Numerical Solution

# cpr.run('height, pressure = ex.squareWave_gmresSolve(domain, p0, pN, n_steps)')

# ncalls  tottime  percall  cumtime  percall filename:lineno(function)
# ---------- N = 101 ----------------
# 8441    0.716    0.000    0.721    0.000 squareWave.py:236(_matvec) <--  class-scope
# 8441    0.644    0.000    0.648    0.000 squareWave.py:280(_matvec) <--  fun-scope 
# ---------- N = 501 ----------------
# 177542   77.865    0.000   77.945    0.000 squareWave.py:236(_matvec) <--  class-scope
# 177542   69.212    0.000   69.293    0.000 squareWave.py:280(_matvec) <--  fun-scope




# height, pressure = ex.squareWave_gmresSolve(domain, p0, pN, n_steps)
# height, pressure = ex.squareWave_pySolve(domain, p0, pN, n_steps)
# TODO: make a wrapper for fd solve
# pressure_fd = ry.solve(domain, height, p0, pN)

# num_err = np.max(np.abs(num_pressure - pressure.ps))
# print("Analytic to Numerical Error: %.3f"%num_err)

# nump_h_title = "Reynolds Numerical Pressure for %s"%height.h_str
# graph.plot_2D_twin(pressure.ps, height.hs, domain.xs, nump_h_title)


#------------------------------------------------------------------------------
# Error

# err = la.norm(pressure_lu.ps - pressure_gmres.ps)
# print("analytic to numerical error %.5f"%(err))



#-------------------------
# Velocity


ys = np.linspace(0, 1.2*height.h_max, domain.Nx)

v_x = 1/(2*domain.eta) * pressure.pxs * (ys**2 - ys*height.hs) + domain.U * (1- ys/height.hs)

v_y = -1/(2*domain.eta) * ( pressure.pxxs * (1/3 * ys**3 - 1/2 * height.hs * ys**2) - 1/2 * domain.U * height.hxs * pressure.pxs * ys**2 - domain.U**2 * domain.eta / height.hs**2 * ys**2 * height.hxs)


phv_title = "Pressure and Velocity for %s"%height.h_str
phv_fun_labels = ['velocity $(v_x, v_y)$', 'pressure $p(x)$', 'height $h(x)$']
phv_ax_labels =  ['x', 'y']

graph.plot_phv(pressure.ps, height.hs, v_x, v_y, domain.xs, ys, phv_title, phv_fun_labels,  phv_ax_labels)

    

#---------------------------------------------------------------------------
# Plotting 


p_h_title = "Pressure (%s) and Height for %s"%(pressure.p_str, height.h_str)
p_h_labels = ["Pressure $p(x)$", "Height $h(x)$", "$x$"]
graph.plot_2D_twin(pressure.ps, height.hs, domain.xs, p_h_title, p_h_labels)

# graph.plot_2D_twin(pressure.ps, height.hs, domain.xs, p_h_title, p_h_labels)

# graph.plot_2D_twin(pressure.ps - pressure_gmres.ps, height.hs, domain.xs, "LU - Gmres Pressure", p_h_labels)




