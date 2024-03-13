#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 22 15:05:40 2023
@author: sarahdennis
"""
import numpy as np

import domain as dm
import pressure_finDiff as fd
import examples as ex
import velocity as vel
import graphics as graph

#------------------------------------------------------------------------------
# Domain & Discretization
#------------------------------------------------------------------------------
# x boundary
x0 = 0
xf = 1
BC = "fixed" # Boundary Condition in x (alt. "periodic")

# surface velocity 
U = 1   #V_x(x,0) = U

# kinematic viscosity
visc = 1   

# Grid size
Nx = 500

domain = dm.Domain(x0, xf, U, Nx, BC)

#------------------------------------------------------------------------------
# Pressure & Height 
#------------------------------------------------------------------------------
# Pressure boundary
p0 = 0

pN = 1

# Height params (see Examples for more)
n = 1

r = 0.25
h_avg = 0.75

x_step = 0.4

h_min = .1
h_max = h_min + 0.5

#------------------------------------------------------------------------------
# Analytic Solutions
#------------------------------------------------------------------------------
height, pressure = ex.flat(domain, p0, pN, h_avg)

# height, pressure = ex.wedge(domain, p0, pN, h_max, h_min)

# height, pressure = ex.corrugated(domain, p0, pN)

# height, pressure = ex.step(domain, p0, pN, x_step, h_avg-r, h_avg+r)

# height, pressure = ex.sawtooth(domain, p0, pN, h_min, h_max, n)
# height, pressure = ex.sawtooth_gk(domain, p0, pN, h_min, h_max)

#------------------------------------------------------------------------------
# Numerical Solutions
#------------------------------------------------------------------------------

# height, pressure = ex.fdSolve(domain, p0, pN)

# Sawtooth

# height, pressure = ex.sawtooth_finDiff(domain, p0, pN, h_min, h_max, n)

# Square wave

# height, pressure  = ex.squareWave_pySolve(domain, p0, pN, n)

# height, pressure = ex.squareWave_schurInvSolve(domain, p0, pN, n, r, h_avg)

# height, pressure = ex.squareWave_schurGmresSolve(domain, p0, pN, n, r, h_avg)

# height, pressure = ex.squareWave_schurLUSolve(domain, p0, pN, n_steps, r, h_avg)

# height, pressure = ex.randSteps_schurLUSolve(domain, p0, pN, h_min, h_max, n)


#------------------------------------------------------------------------------
# Pressure & Height plotting 
#------------------------------------------------------------------------------
ph_title = "Pressure & Height: \n %s"%(pressure.p_str)
ph_labels = ["Pressure $p(x)$", "Height $h(x)$", "$x$"]
graph.plot_2D_twin(pressure.ps, height.hs, domain.xs, ph_title, ph_labels)


#------------------------------------------------------------------------------
# Velocity
#------------------------------------------------------------------------------
velocity = vel.Velocity(domain, height, pressure)

v_title = "Velocity: \n %s"%(pressure.p_str)
v_labels = ['velocity $(v_x, v_y)$', 'pressure $p(x)$', 'height $h(x)$']
v_ax_labels =  ['$x$', '$y$']

graph.plot_stream(velocity.vx, velocity.vy, domain.xs, domain.ys, v_title, v_ax_labels)

#------------------------------------------------------------------------------
# Error
#------------------------------------------------------------------------------
# num_anl_title = "Numerical and Analytic Pressure for %s"%height.h_str
# num_anl_labels = [pressure.p_str, "finDiff"]
# num_anl_axis = ["$x$", "Error"]
# graph.plot_2D_multi([pressure.ps, num_pressure], domain.xs, num_anl_title, num_anl_labels, num_anl_axis)

# max_num_err = np.max(np.abs(num_pressure - pressure.ps))
# print("Analytic to Numerical Error: %.3f"%max_num_err)

# num_err_title = "Numerical Error for %s"%height.h_str
# num_err_axis = ["$x$", "Pressure $p$"]
# graph.plot_2D(pressure.ps - num_pressure, domain.xs, num_err_title, num_err_axis)
