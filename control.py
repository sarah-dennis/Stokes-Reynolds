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
# x boundary
x0 = 0     
xf = 1      

BC = "fixed" # Boundary Condition in x (alt. "periodic")

# surface velocity 
U = 0     #V_x(x,0) = U
    
# viscosity
visc = 1     

# Grid size
#TODO: Nx used in discretization & plotting. Using Ny = Nx.
Nx = 1000 

domain = dm.Domain(x0, xf, visc, U, Nx, BC)

#------------------------------------------------------------------------------
# Height & Pressure
#-------------------------------------------------------------------------------
p0 = 2 #
pN = 1

n_steps = 5

r=0.02

h_avg=0.1

#---------------------------------------------------------------------------
# Analytic Solutions

# height, pressure, velocity = ex.flat(domain, p0, pN, h_avg)
# height, pressure, velocity = ex.wedge(domain, p0, pN)
# height, pressure, velocity = ex.corrugated(domain, p0, pN)
height, pressure, velocity = ex.step(domain, p0, pN, r, h_avg)


# height, pressure, velocity = ex.squareWave_schurInvSolve(domain, p0, pN, n_steps, r, h_avg)
height, pressure, velocity = ex.squareWave_schurLUSolve(domain, p0, pN, n_steps, r, h_avg)

# height, pressure, velocity = ex.randRectWave_schurLUSolve(domain, p0, pN, n_steps, r, h_avg)

#---------------------------------------------------------------------------
# Numerical Solution


# cpr.run('height, pressure = ex.squareWave_gmresSolve(domain, p0, pN, n_steps)')

# Square wave: numerical solves
# height, pressure, velocity = ex.squareWave_gmresSolve(domain, p0, pN, n_steps)
# height, pressure, velocity = ex.squareWave_pySolve(domain, p0, pN, n_steps)

# Random height: finite difference sovle
# height, pressure, velocity = ex.randGrid(domain, p0, pN)

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
        
phv_title = "Pressure and Velocity for %s"%height.h_str
phv_fun_labels = ['velocity $(v_x, v_y)$', 'pressure $p(x)$', 'height $h(x)$']
phv_ax_labels =  ['$x$', '$y$']

graph.plot_phv(pressure.ps, height.hs, velocity.vx, velocity.vy, domain.xs, domain.ys, phv_title, phv_fun_labels,  phv_ax_labels)


