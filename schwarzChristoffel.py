#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 30 11:05:44 2023

@author: sarahdennis
"""
import numpy as np
import graphics
import cmath

 
#------------------------------------------------------------------------------
# Pentagon
#------------------------------------------------------------------------------
# -  M ----------- L -
# h1 |             | |
# -  N ---- P      | h2
#           |      | |
#           R ---- S -     
#    |--l1--|--l2--|  

# flow:  <-#-U-- 
U=1

l1 = 4 
l2 = 6 

h1 = 1 
h2 = 2 

# h2 > h1, a = (H/h)**2
#------------------------------------------------------------------------------    
# t : Upper half plane
#             .
#             |
#             |
#             |
# -infty ---- 0 ---- 1 --- a ---- +infty
#  <-- L --- N=M --- P --- R ---- S -->

# t = x + iy

Nx = 100 # f(t_Nx) ~ z_Ny
Ny = 300 # f(t_Ny) ~ z_Nx

t_xMin = -100
t_xMax = 100
t_yMin = 0
t_yMax = 200

t_dx = (t_xMax - t_xMin)/(Nx-1)
t_dy = (t_yMax - t_yMin)/(Ny-1)

t_xs = np.zeros(Nx)
t_ys = np.zeros(Ny)

eps_x = 1/t_dx #avoid t=0
eps_y = 1/t_dy

for i in range (Nx):
    t_x = t_xMin + i*t_dx
    if t_x == 0:
        t_xs[i]= t_x+t_dx*eps_x
    else:     
        t_xs[i] = t_x
for j in range (Ny):
    t_y = t_yMin + j*t_dy
    if t_y == 0:
        t_ys[j] = t_y + t_dy*eps_y
    else:  
        t_ys[j] = t_y
     
#------------------------------------------------------------------------------   
# # Schwarz-Christoffel Mapping z = f(t)
#------------------------------------------------------------------------------       
a = (h2/h1)**2
b1 = h2/np.pi 
b2 = h1/np.pi
K = np.pi/(U*h2)
S = (U*h2)/np.pi

def W(t): #W = phi + i psi
    return S*cmath.log(t)

def f(w):
    
    eKw = np.exp(K*w.real)
    sinKw = cmath.sin(K*w.imag)
    cosKw = cmath.cos(K*w.imag)

    u1 = ( 2 * eKw   * cosKw - a-1)/(a-1)
    v1 = ( 2 * eKw   * sinKw      )/(a-1)
    u2 = (-2 * a/eKw * cosKw + a+1)/(a-1)
    v2 = ( 2 * a/eKw * sinKw      )/(a-1)
    
    alph1 = cmath.sqrt((1+u1)**2 + v1**2)
    beta1 = cmath.sqrt((1-u1)**2 + v1**2)
    alph2 = cmath.sqrt((1+u2)**2 + v2**2)
    beta2 = cmath.sqrt((1-u2)**2 + v2**2)

    acosh1_real = cmath.acosh((alph1 + beta1)/2)
    acosh1_imag = cmath.acosh((alph1 - beta1)/2) 
    acosh2_real = cmath.acosh((alph2 + beta2)/2)
    acosh2_imag = cmath.acosh((alph2 - beta2)/2)

    f_real = b1 * acosh1_real - b2 * acosh2_real
    f_imag = b1 * acosh1_imag - b2 * acosh2_imag - complex(0,(h2-h1))
    return f_real + f_imag


f_zs = np.zeros(Nx * Ny, complex)
for i in range(Nx):
    for j in range(Ny):
        t = complex(t_xs[i], t_ys[j])
        f_zs[i*Ny + j] = f(W(t))

#------------------------------------------------------------------------------
# Stream Function  s := phi(z)
#------------------------------------------------------------------------------
def phi(z): #stream
    return (S/2j)*cmath.log(z/z.conjugate())  
 
stream_zs = np.zeros(Nx * Ny, complex)
for k in range(Nx * Ny):
    stream_zs[k] = phi(f_zs[k])
    
#------------------------------------------------------------------------------
# Velocity Potential Function p := psi(z)
#------------------------------------------------------------------------------
def psi(z):
    return (S/2)*cmath.log(z*z.conjugate())

velpot_zs = np.zeros(Nx * Ny, complex)
for k in range(Nx * Ny):
    velpot_zs[k] = psi(f_zs[k])

print("Complete: starting plot")
#------------------------------------------------------------------------------
# Plotting
#------------------------------------------------------------------------------
stream_title = "Stream $t \in [%d, %d]x[0, %d]$"%(t_xMin, t_xMax, t_yMax)
velpot_title = "Velocity Potential $t \in [%d, %d]x[0, %d]$"%(t_xMin, t_xMax, t_yMax)
ax_labels = ["x", "y"]

f_xs = f_zs.real
f_ys = f_zs.imag
stream_xs = stream_zs.real
stream_ys = stream_zs.imag
velpot_xs = velpot_zs.real
velpot_ys = velpot_zs.imag 

zip_fxy = np.array(sorted(zip(f_xs, f_ys, stream_xs, stream_ys, velpot_xs, velpot_ys))).T  #sort in x
f_xs = zip_fxy[0]
f_ys = zip_fxy[1]
stream_xs = zip_fxy[2]
stream_ys = zip_fxy[3]
velpot_xs = zip_fxy[4]
velpot_ys = zip_fxy[5]

graphics.plot_quivers_flat(stream_xs, stream_ys, f_xs, f_ys, stream_title, ax_labels)

# graphics.plot_quivers_flat(velpot_xs, velpot_ys, f_xs, f_ys, velpot_title, ax_labels)