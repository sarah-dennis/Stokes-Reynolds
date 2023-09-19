#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 22 15:05:40 2023

@author: sarahdennis
"""

import numpy as np

import cProfile as cpr

import domain as dm
import finDiff_1D as fd
import examples_1D as ex
import _graphics as graph


#------------------------------------------------------------------------------
# Domain & Discretization
#------------------------------------------------------------------------------
x0 = 0      # left boundary 
xf = 1      # right boundary

U = 5     # surface velocity

eta = 1     # viscosity


Nx = 1000   # Number of Grid points

BC = "fixed" # Boundary Condition in x (alt. "periodic")

domain = dm.Domain(x0, xf, eta, U, Nx, BC)

#------------------------------------------------------------------------------
# Height & Pressure
#-------------------------------------------------------------------------------
p0 = 0
pN = 0

n_steps = 3
r=0.028
h_avg=0.1

#---------------------------------------------------------------------------
# Analytic Solution

# height, pressure = ex.flat(domain, p0, pN)
height, pressure = ex.wedge(domain, p0, pN)
# height, pressure = ex.corrugated(domain, p0, pN)
# height, pressure = ex.step(domain, p0, pN)


# height, pressure = ex.squareWave_schurInvSolve(domain, p0, pN, n_steps, r, h_avg)
height, pressure = ex.squareWave_schurLUSolve(domain, p0, pN, n_steps, r, h_avg)



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



# Square wave: numerical matrix solves
# height, pressure = ex.squareWave_gmresSolve(domain, p0, pN, n_steps)
# height, pressure = ex.squareWave_pySolve(domain, p0, pN, n_steps)

# Random height: finite difference sovle
# height, pressure = ex.randGrid(domain, p0, pN)


#---------------------------------------------------------------------------
# Plotting 

p_h_title = "Pressure (%s) and Height for %s"%(pressure.p_str, height.h_str)
p_h_labels = ["Pressure $p(x)$", "Height $h(x)$", "$x$"]
graph.plot_2D_twin(pressure.ps, height.hs, domain.xs, p_h_title, p_h_labels)

#------------------------------------------------------------------------------
# Error
# num_anl_title = "Numerical and Analytic Pressure for %s"%height.h_str
# num_anl_labels = [pressure.p_str, "finDiff"]
# num_anl_axis = ["$x$", "Error"]
# graph.plot_2D_multi([pressure.ps, num_pressure], domain.xs, num_anl_title, num_anl_labels, num_anl_axis)



# max_num_err = np.max(np.abs(num_pressure - pressure.ps))
# print("Analytic to Numerical Error: %.3f"%max_num_err)
# num_err_title = "Numerical Error for %s"%height.h_str
# num_err_axis = ["$x$", "Pressure $p$"]

# graph.plot_2D(pressure.ps - num_pressure, domain.xs, num_err_title, num_err_axis)

#-------------------------
# Velocity

ys = np.linspace(0, 1.2*height.h_max, domain.Nx)

vx = np.zeros((domain.Nx, domain.Nx))

for j in range(domain.Nx):
    for i in range(domain.Nx):
        vx[j,i] = 1/(2*domain.eta) * pressure.pxs[i] * (ys[j]**2 - height.hs[i]*ys[j]) + domain.U * (1 - ys[j]/height.hs[i])

phv_title = "Pressure and Velocity for %s"%height.h_str
phv_fun_labels = ['velocity $(v_x, v_y)$', 'pressure $p(x)$', 'height $h(x)$']
phv_ax_labels =  ['x', 'y']

graph.plot_phv(pressure.ps, height.hs, vx, vy, domain.xs, ys, phv_title, phv_fun_labels,  phv_ax_labels)

    




