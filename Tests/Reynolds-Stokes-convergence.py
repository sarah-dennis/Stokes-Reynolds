#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr  1 12:28:25 2024

@author: sarahdennis
"""

import numpy as np
import graphics as graph

import domain as dm
import heights as hgt

import pressures as reyns
import velocity as vel

import Stokes as stokes

#------------------------------------------------------------------------------
# Domain 
#------------------------------------------------------------------------------
# x boundary
x0 = 0
xN = 1

# pressure boundary
BC = "fixed" 
p0 = 0
pN = 0

# surface velocity 
U = 1   #V_x(x,0) = U

# kinematic viscosity
visc = 1   

# Numerical grid size
N = 1000

Nx = N+1
Ny = (2*N)+1

domain = dm.Domain(x0, xN, U, Nx, Ny, BC)

#------------------------------------------------------------------------------
# Height
#------------------------------------------------------------------------------
n = 2 # number of linear components > 0

h_min = 1e-3 *domain.dx
h_max = 2

xi = domain.x0 + np.arange(0, n+1) * (domain.xf - domain.x0)/n

hi = np.array([h_min, h_max, h_min])

height = hgt.SawtoothHeight(domain, xi, hi)

#------------------------------------------------------------------------------
# Reynolds solution
#------------------------------------------------------------------------------
reyn_pressure = reyns.SawtoothPressure(domain, height, p0, pN)

reyn_vel = vel.Velocity(domain, height, reyn_pressure)

reyn_u, reyn_v = np.flip(reyn_vel.vx), np.flip(reyn_vel.vy)
reyn_uv = np.dstack((reyn_u, reyn_v))

hs = height.h_max - height.hs


reyn_title = "Velocity: \n %s"%(reyn_pressure.p_str)
ax_labels_vel =  ['$x$', '$y$']

graph.plot_stream_height(reyn_u, reyn_v, hs, domain.xs, domain.ys, reyn_title, ax_labels_vel)


#------------------------------------------------------------------------------
# Stokes solution
#------------------------------------------------------------------------------

filename = "stokes_N%d_spLU.csv" % N

stokes_u_flat, stokes_v_flat, psi, past_iters = stokes.read_solution(filename, Nx*Ny)

stokes_stream = psi.reshape(Ny, Nx)
stokes_u = stokes_u_flat.reshape(Ny, Nx)
stokes_v = stokes_v_flat.reshape(Ny, Nx)
stokes_uv = np.dstack((stokes_u, stokes_v))

stokes_title_stream = "Stream: \n Biharmonic Stokes $k=%d$"%past_iters
ax_labels_stream = ['$\psi$', '$x$', '$y$']
# graph.plot_heat_contour(stokes_stream, domain.xs, domain.ys, stokes_title_stream, ax_labels_stream)

stokes_title_vel = "Velocity: \n Biharmonic Stokes $k=%d$"%past_iters
graph.plot_stream_height(stokes_u, stokes_v, hs, domain.xs, domain.ys, stokes_title_vel, ax_labels_vel)


#------------------------------------------------------------------------------
# Error calculations
#------------------------------------------------------------------------------

def l1pctError(f_u, f_v, g_u, g_v, Nx, Ny):
    max_err = 0.0
    norm_err = np.zeros((Ny, Nx))

    for j in range(Ny):
        for i in range(Nx):
                        
            fu = f_u[j,i]
            
            fv = f_v[j,i]
            
            gu = g_u[j,i]
            
            gv = g_v[j,i]
            
            err_ij = abs(fu - gu) + abs(fv - gv)
            
            norm_ij = abs(gu) + abs(gv)
            
            if norm_ij == 0:
                norm_err_ij = 0
            else:
                norm_err_ij = 100*err_ij/norm_ij
            norm_err[j,i] = norm_err_ij
            max_err = max(norm_err_ij, max_err)
    return norm_err, max_err

def l2pctError(f_u, f_v, g_u, g_v, Nx, Ny):
    max_err = 0.0
    
    norm_err = np.zeros((Ny, Nx))
    for j in range(Ny):
        for i in range(Nx):
            fu = f_u[j,i]
            
            fv = f_v[j,i]
            
            gu = g_u[j,i]
            
            gv = g_v[j,i]
           
            err_ij = np.sqrt((fu - gu)**2 + (fv - gv)**2)

            norm_ij = np.sqrt(gu**2 + gv**2)
            
            if norm_ij == 0:
                norm_err_ij = 0
            else:
                norm_err_ij = 100* err_ij/norm_ij
          
            norm_err[j,i] = norm_err_ij
            
            max_err = max(norm_err_ij, max_err)
    return norm_err, max_err


def lInfpctError(f_u, f_v, g_u, g_v, Nx, Ny):

    max_err = 0.0
    norm_err = np.zeros((Ny, Nx))

    for j in range(Ny):
        for i in range(Nx):
            fu = f_u[j,i]
            
            fv = f_v[j,i]
            
            gu = g_u[j,i]
            
            gv = g_v[j,i]
           
            err_ij = max(abs(fu - gu), abs(fv - gv))
            
            norm_ij = max(abs(gu), abs(gv))
            
            if norm_ij == 0:
                norm_err_ij = 0
            else:
                norm_err_ij = 100* err_ij/norm_ij
            norm_err[j,i] = norm_err_ij
            max_err = max(norm_err_ij, max_err)
    return norm_err, max_err


l1_err, l1_maxerr = l1pctError(reyn_u, reyn_v, stokes_u, stokes_v, Nx, Ny)
l2_err, l2_maxerr = l2pctError(reyn_u, reyn_v, stokes_u, stokes_v, Nx, Ny)
linf_err, linf_maxerr= lInfpctError(reyn_u, reyn_v, stokes_u, stokes_v, Nx, Ny)

title_l1_err = "$l_1$ Percent error (Reynolds : Stokes)"
ax_labels_l1_err = ['$|u_r - u_s|/|u_s|$','x','y']
graph.plot_heat_contour(l2_err, domain.xs, domain.ys, title_l1_err, ax_labels_l1_err)


title_l2_err = "$l_2$ Percent error"
ax_labels_l2_err = ['$|u_r - u_s|_2/|u_s|_2$','x','y']
graph.plot_heat_contour(l2_err, domain.xs, domain.ys, title_l2_err, ax_labels_l2_err)


title_linf_err = "$l_\infty$ Percent error"
ax_labels_linf_err = ['$|u_r - u_s|_\infty/|u_s|_\infty$','x','y']
graph.plot_heat_contour(l2_err, domain.xs, domain.ys, title_linf_err, ax_labels_linf_err)
