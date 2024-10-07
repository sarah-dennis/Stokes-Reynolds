#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr  1 12:28:25 2024

@author: sarahdennis
"""

import numpy as np
import graphics as graph

import domain as dm
import reyn_heights

import reyn_examples as rex
import reyn_velocity as vel

import stokes_control as stokes
import stokes_readwrite as rw

import stokes_examples as sex
import error 
#------------------------------------------------------------------------------
# Domain 
#------------------------------------------------------------------------------
# x boundary
x0 = 0
xf = 1

# pressure boundary
BC = "fixed" 
p0 = 0
pN = 0

# surface velocity 
U = 1   #V_x(x,0) = U

# kinematic viscosity
visc = 1   

# Numerical grid size
N = 200

Nx = N+1
Ny = (2*N)+1

#------------------------------------------------------------------------------
# Reynolds solution
#------------------------------------------------------------------------------
# reyn_pressure = rex.PWA_Sawtooth(height, p0, pN)
reyn_pressure = rex.StepWave_Ex3()
height = reyn_pressure.height
reyn_vel = vel.Velocity(height, reyn_pressure.ps)

reyn_u, reyn_v = np.flip(reyn_vel.vx, axis=0), np.flip(reyn_vel.vy, axis=0)
reyn_uv = np.dstack((reyn_u, reyn_v))

hs = height.h_max - height.hs


reyn_title = "Velocity: \n Reynolds Stepwave"
ax_labels_vel =  ['$x$', '$y$']

graph.plot_stream_height(reyn_u, reyn_v, hs, height.xs, height.ys, reyn_title, ax_labels_vel)


#------------------------------------------------------------------------------
# Stokes solution
#------------------------------------------------------------------------------

# tri = sex.biswasEx(N)

# stokes_u_flat, stokes_v_flat, psi, past_iters = rw.read_solution(tri.filename+".csv", Nx*Ny)

# stokes_stream = psi.reshape(Ny, Nx)
# stokes_u = stokes_u_flat.reshape(Ny, Nx)
# stokes_v = stokes_v_flat.reshape(Ny, Nx)
# stokes_uv = np.dstack((stokes_u, stokes_v))

# stokes_title_vel = "Velocity: \n Stokes $N=%d$"%N
# graph.plot_stream_height(stokes_u, stokes_v, hs, tri.xs, tri.ys, stokes_title_vel, ax_labels_vel)


#------------------------------------------------------------------------------
# Error calculations
#------------------------------------------------------------------------------

# l1_err, l1_maxerr = error.l1pctError(reyn_u, reyn_v, stokes_u, stokes_v, Nx, Ny)
# l2_err, l2_maxerr = error.l2pctError(reyn_u, reyn_v, stokes_u, stokes_v, Nx, Ny)
# linf_err, linf_maxerr= error.lInfpctError(reyn_u, reyn_v, stokes_u, stokes_v, Nx, Ny)

# print(l1_maxerr)
# print(l2_maxerr)
# print(linf_maxerr)
