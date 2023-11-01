#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 22 15:05:40 2023
@author: sarahdennis
"""

import numpy as np

import domain as dm
import finDiff_1D as fd
import examples_1D as ex
import velocity as vel
import _graphics as graph

#------------------------------------------------------------------------------
# Domain & Discretization
#------------------------------------------------------------------------------
# x boundary
x0 = 0
xf = 5    
BC = "fixed" # Boundary Condition in x (alt. "periodic")

# surface velocity 
U = 2   #V_x(x,0) = U

# kinematic viscosity
visc = 1   

# Grid size
Nx = 500

domain = dm.Domain(x0, xf, visc, U, Nx, BC)

#------------------------------------------------------------------------------
# Pressure & Height 
#------------------------------------------------------------------------------
# Pressure boundary
p0 = 0
pN = 0

# Height params (see Examples for more)
n_steps = 57
r = 0.25
h_avg = 0.75

x_step = 1

#------------------------------------------------------------------------------
# Analytic Solutions
#------------------------------------------------------------------------------
# height, pressure = ex.flat(domain, p0, pN, h_avg)
# height, pressure = ex.wedge(domain, p0, pN)
height, pressure = ex.corrugated(domain, p0, pN)

# height, pressure = ex.step(domain, p0, pN, x_step, r, h_avg)

#------------------------------------------------------------------------------
# Numerical Solutions
#------------------------------------------------------------------------------
# Square wave: numerical solves
# height, pressure = ex.squareWave_schurLUSolve(domain, p0, pN, n_steps, r, h_avg)

# height, pressure = ex.squareWave_schurInvSolve(domain, p0, pN, n_steps, r, h_avg)

# height, pressure = ex.squareWave_gmresSolve(domain, p0, pN, n_steps, r, h_avg)
# height, pressure  = ex.squareWave_pySolve(domain, p0, pN, n_steps)
# h_steps = [2, 2, 1, 1]
# height, pressure = ex.mySteps_schurLUSolve(domain, p0, pN, h_steps)


#------------------------------------------------------------------------------
# Pressure & Height plotting 
#------------------------------------------------------------------------------
p_h_title = "Pressure (%s) and Height for %s"%(pressure.p_str, height.h_str)
p_h_labels = ["Pressure $p(x)$", "Height $h(x)$", "$x$"]
graph.plot_2D_twin(pressure.ps, height.hs, domain.xs, p_h_title, p_h_labels)


#------------------------------------------------------------------------------
# Velocity
#------------------------------------------------------------------------------
velocity = vel.Velocity(domain, height, pressure)

phv_title = "Pressure and Velocity for %s"%height.h_str
phv_fun_labels = ['velocity $(v_x, v_y)$', 'pressure $p(x)$', 'height $h(x)$']
phv_ax_labels =  ['$x$', '$y$']

graph.plot_phv(pressure.ps, height.hs, velocity.vx, velocity.vy, domain.xs, domain.ys, phv_title, phv_fun_labels,  phv_ax_labels)

# Inlet velocity
velocity.plot_vx_x0(domain, 0)
# velocity.plot_vy_x0(domain, 0)

#Step profile 
# index = domain.get_index(x_step)
# velocity.plot_vx_x0(domain, index)
# velocity.plot_vy_x0(domain, index)

# Outlet velocity
velocity.plot_vx_x0(domain, -1)
# velocity.plot_vy_x0(domain, -1)



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
