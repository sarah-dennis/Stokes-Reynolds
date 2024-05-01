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
N = 240

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
hs = height.h_max - height.hs


reyn_title = "Velocity: \n %s"%(reyn_pressure.p_str)
ax_labels =  ['$x$', '$y$']

graph.plot_stream_height(reyn_u, reyn_v, hs, domain.xs, domain.ys, reyn_title, ax_labels)


#------------------------------------------------------------------------------
# Stokes solution
#------------------------------------------------------------------------------

filename = "biharmonicStokes_N%d.csv" % N

stokes_u_flat, stokes_v_flat, psi, past_iters = stokes.read_solution(filename, Nx*Ny)

stokes_stream = psi.reshape(Ny,Nx)
stokes_u = stokes_u_flat.reshape(Ny, Nx)
stokes_v = stokes_v_flat.reshape(Ny, Nx)

stokes_title_stream = "Stream: \n Biharmonic Stokes $k=%d$"%past_iters
ax_labels_stream = ['$\psi$', '$x$', '$y$']
graph.plot_heat_contour(stokes_stream, domain.xs, domain.ys, stokes_title_stream, ax_labels_stream)

stokes_title_vel = "Velocity: \n Biharmonic Stokes $k=%d$"%past_iters
graph.plot_stream_height(stokes_u, stokes_v, hs, domain.xs, domain.ys, stokes_title_vel, ax_labels)


#------------------------------------------------------------------------------
# Error calculations
#------------------------------------------------------------------------------
def L1Max(us, vs, Nx, Ny):
    max_n = 0
    for j in range(Ny):
        for i in range(Nx):
            u = abs(us[j,i])
            v = abs(vs[j,i])
            n = u + v
            max_n = max(n, max_n)
    return max_n

def LinfMax(us, vs, Nx, Ny):
    max_n = 0
    for j in range(Ny):
        for i in range(Nx):
            u = abs(us[j,i])
            v = abs(vs[j,i])
            n = max(u,v)
            max_n = max(n, max_n)
    return max_n


def L2Max(us, vs, Nx, Ny):
    max_n = 0
    for j in range(Ny):
        for i in range(Nx):
            u = abs(us[j,i])
            v = abs(vs[j,i])
            n = (u**2 + v**2)**(1/2)
            max_n = max(n, max_n)
    return max_n


def L1Error(us_reyn, vs_reyn, us_stokes, vs_stokes, Nx, Ny):
    max_err = 0.0
    
    for j in range(Ny):
        for i in range(Nx):
                        
            u_r = us_reyn[j,i]
            
            u_s = us_stokes[j,i]
            
            v_r = vs_reyn[j,i]
            
            v_s = vs_stokes[j,i]
            
            err_ij = abs(u_r - u_s) + abs(v_r - v_s)
            
            # norm_ij = abs(u_s) + abs(v_s)
            
            # if norm_ij == 0:
            #     continue
            # else:
            #     normErr_ij = err_ij/norm_ij
            #     max_err = max(normErr_ij, max_err)
            max_err = max(err_ij, max_err)
    return max_err

def L2Error(us_reyn, vs_reyn, vs_stokes, us_stokes, Nx, Ny):
    max_err = 0.0
    
    for j in range(Ny):
        for i in range(Nx):
                        
            u_r = us_reyn[j,i]
            
            u_s = us_stokes[j,i]
            
            v_r = vs_reyn[j,i]
            
            v_s = vs_stokes[j,i]
            
            err_ij = (abs(u_r - u_s)**2 + abs(v_r - v_s)**2)**(1/2)
            
            # norm_ij = (u_s**2 + v_s**2)**(1/2)
            
            # normErr_ij = err_ij/norm_ij
            
            # if norm_ij == 0:
            #     continue
            # else:
            #     normErr_ij = err_ij/norm_ij
            #     max_err = max(normErr_ij, max_err)
            max_err = max(err_ij, max_err)
    return max_err


def LinfError(us_reyn, vs_reyn, vs_stokes, us_stokes, Nx, Ny):

    max_err = 0.0
    for j in range(Ny):
        for i in range(Nx):
            u_r = us_reyn[j,i]
            
            u_s = us_stokes[j,i]
            
            v_r = vs_reyn[j,i]
            
            v_s = vs_stokes[j,i]
           
            err_ij = max(abs(u_r - u_s), abs(v_r - v_s))
            
            # norm_ij = max(abs(u_s), abs(v_s))
            
            # normErr_ij = err_ij/norm_ij
            
            # if norm_ij == 0:
            #     continue
            # else:
            #     normErr_ij = err_ij/norm_ij
            #     max_err = max(normErr_ij, max_err)
            max_err = max(err_ij, max_err)
    return max_err

l1max = L1Max(stokes_u, stokes_v, Nx, Ny)
l2max = L2Max(stokes_u, stokes_v, Nx, Ny)
linfmax = LinfMax(stokes_u, stokes_v, Nx, Ny)


Y = cdist(XA, XB, 'minkowski', p=2.)

Y = cdist(XA, XB, 'seuclidean', V=None)


l1_pcterr = 100*L1Error(reyn_u, reyn_v, stokes_u, stokes_v, Nx, Ny)/l1max
l2_pcterr = 100*L2Error(reyn_u, reyn_v, stokes_u, stokes_v, Nx, Ny)/l2max
linf_pcterr = 100*LinfError(reyn_u, reyn_v, stokes_u, stokes_v, Nx, Ny)/linfmax

print('l1 pct error: %.4f'%l1_pcterr)
print('l2 pct error: %.4f'%l2_pcterr)
print('l* pct error: %.4f'%linf_pcterr)
