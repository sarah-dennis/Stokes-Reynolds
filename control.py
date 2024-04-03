#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 22 15:05:40 2023
@author: sarahdennis
"""
import numpy as np
import csv
import domain as dm
import examples as ex
import velocity as vel
import graphics as graph

#------------------------------------------------------------------------------
# Domain & Discretization
#------------------------------------------------------------------------------
# x boundary
x0 = 0
xN = 1
BC = "fixed" # Boundary Condition in x (alt. "periodic")

# surface velocity 
U = 1   #V_x(x,0) = U

# kinematic viscosity
visc = 1   

# Grid size
N = 2000

domain = dm.Domain(x0, xN, U, N, N, BC)

#------------------------------------------------------------------------------
# Pressure & Height 
#------------------------------------------------------------------------------
# Pressure boundary
p0 = 0

pN = 0

# Height params (see Examples for more)
n = 4

x_step = 0.5  #x0 < step < xN

r = 1

h_min = .1
h_max = h_min + r

#------------------------------------------------------------------------------
# Analytic Solutions
#------------------------------------------------------------------------------
# height, al_pressure = ex.constant(domain, p0, pN, h_avg)

# height, al_pressure = ex.corrugated(domain, p0, pN)

# height, al_pressure = ex.step(domain, p0, pN, x_step, h_min, h_max)


# height, al_pressure= ex.sawtooth(domain, p0, pN, h_min, h_max, n)

height, al_pressure = ex.sawtoothRand(domain, p0, pN, h_min, h_max, n)

#------------------------------------------------------------------------------
# Numerical Solutions
#------------------------------------------------------------------------------

fd_pressure = ex.fdSolve(domain, p0, pN, height)

# Square wave

# height, pressure  = ex.squareWave_pySolve(domain, p0, pN, n)

# height, pressure = ex.squareWave_schurInvSolve(domain, p0, pN, n, r, h_avg)

# height, pressure = ex.squareWave_schurGmresSolve(domain, p0, pN, n, r, h_avg)

# height, pressure = ex.squareWave_schurLUSolve(domain, p0, pN, n_steps, r, h_avg)

# height, pressure = ex.randSteps_schurLUSolve(domain, p0, pN, h_min, h_max, n)

#------------------------------------------------------------------------------
# Pressure & Height plotting 
#------------------------------------------------------------------------------
def plot_ph(domain, pressure, height):
    ph_title = "Pressure & Height: \n %s"%(pressure.p_str)
    ph_labels = ["Pressure $p(x)$", "Height $h(x)$", "$x$"]
    graph.plot_2D_twin(pressure.ps, height.hs, domain.xs, ph_title, ph_labels)

plot_ph(domain, fd_pressure, height)
plot_ph(domain, al_pressure, height)
#------------------------------------------------------------------------------
# Velocity
#------------------------------------------------------------------------------

def plot_v(domain, pressure, height):
    velocity = vel.Velocity(domain, height, pressure)

    v_title = "Velocity: \n %s"%(pressure.p_str)
    v_ax_labels =  ['$x$', '$y$']

    graph.plot_stream(velocity.vx, velocity.vy, domain.xs, domain.ys, v_title, v_ax_labels)

plot_v(domain, fd_pressure, height)
plot_v(domain, al_pressure, height)

#------------------------------------------------------------------------------
# Error
#------------------------------------------------------------------------------

infNorm_err = np.max(np.abs(al_pressure.ps - fd_pressure.ps))
print("Analytic to Numerical Error: %.8f"%infNorm_err)

num_err_title = "Numerical Error for %s"%height.h_str
num_err_axis = ["$x$", "Error (Analytic - Numerical)"]
graph.plot_2D(al_pressure.ps - fd_pressure.ps, domain.xs, num_err_title, num_err_axis)

#------------------------------------------------------------------------------
# CSV writing
#------------------------------------------------------------------------------

def write_solution(filename, length, fd_ps, al_ps, N, err):

    with open(filename, 'w', newline='') as file:
        writer = csv.writer(file, delimiter=' ')
        writer.writerow(["N=%d"%N])
        writer.writerow(["err:", err])
        writer.writerow(["finDiff", "analytic"])
        
        for i in range(length):
            writer.writerow([fd_ps[i], al_ps[i]])
        
        print("  saved csv")
        file.close()


# filename = "%d_linear_pressure_%d"%(n, N)

# write_solution(filename, domain.Nx, fd_pressure.ps, al_pressure.ps, N, infNorm_err)

