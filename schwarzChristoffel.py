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
# t : Upper half plane
#             .
#             |
#             |
#             |
# -infty ---- 0 ---- 1 --- a ---- +infty
#  <-- L --- N=M --- P --- R ---- S -->

# t = x + iy

Nx = 300
Ny = 100

t_xMin = -500
t_xMax = 500
t_yMin = 0
t_yMax = 100

t_dx = (t_xMax - t_xMin)/(Nx-1)
t_dy = (t_yMax - t_yMin)/(Ny-1)

t_xs = np.zeros(Nx)
t_ys = np.zeros(Ny)

for i in range (Nx):
    t_x = t_xMin + i*t_dx
    if t_x == 0:
        t_xs[i]= t_xMin + i*(t_dx/2)
    else:     
        t_xs[i] = t_x
for j in range (Ny):
    t_y = t_yMin + j*t_dy
    if t_y == 0:
        t_ys[j] = t_yMin + j*(t_dy/2)
    else:  
        t_ys[j] = t_y
    
ts = np.zeros(Nx * Ny, complex)
for i in range (Nx):
    for j in range (Ny):
        ts[i*Ny+j] = complex(t_xs[i], t_ys[j])
            
 # ts[k] = complex(t_xs[k//t_Ny], t_ys[k%t_Ny]) 
 
 
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
U=0.5

l1 = 4 # dist N -> P
l2 = 6 # dist R -> S

h1 = 1 #dist N -> M
h2 = 2 #dist S -> L

# Step height : 
#    h2 > h1, a = (H/h)**2
     
#------------------------------------------------------------------------------   
# # Schwarz-Christoffel Mapping z = f(t)
#------------------------------------------------------------------------------       
a = (h2/h1)**2
b1 = h2/np.pi 
b2 = h1/np.pi


S = (U*h2)/np.pi
K = np.pi/(U*h2)

def W(t): #W = phi + i psi
    return S*cmath.log(t)

def f(w):
    psi_t = w.real
    phi_t = w.imag

    u1 = (2*np.exp(K*psi_t) * cmath.cos(K*phi_t) - (a+1)) /(a-1)
    v1 = (2*np.exp(K*psi_t) * cmath.sin(K*phi_t)) /(a-1)
    
    alph1 = cmath.sqrt((1+u1)**2 + v1**2)
    beta1 = cmath.sqrt((1-u1)**2 + v1**2)

    acosh1_real = cmath.acosh(0.5*(alph1 + beta1))
    acosh1_imag = cmath.acosh(0.5*(alph1 - beta1)) 
    
    u2 = (-2*a*np.exp(-K*psi_t) * cmath.cos(K*phi_t) + (a+1)) /(a-1)
    v2 = (2*a*np.exp(-K*psi_t) * cmath.sin(K*phi_t)) /(a-1)

    alph2 = cmath.sqrt((1+u2)**2 + v2**2)
    beta2 = cmath.sqrt((1-u2)**2 + v2**2)
    
    acosh2_real = cmath.acosh(0.5*(alph2 + beta2))
    acosh2_imag = cmath.acosh(0.5*(alph2 - beta2))

    f_real = b1 * acosh1_real - b2 * acosh2_real
    f_imag = b1 * acosh1_imag - b2 * acosh2_imag - complex(0,(h2-h1))
    
    return complex(f_real.real, f_imag.imag)


f_zs = np.zeros(Nx * Ny, complex)

for k in range(Nx*Ny):
    f_zs[k] = f(W(ts[k]))
    
f_xs = f_zs.real
f_ys= f_zs.imag

#------------------------------------------------------------------------------
# Stream Function Mapping   s := phi(z)
#------------------------------------------------------------------------------
def stream(z): #stream
     return (S/2j)*cmath.log(z/z.conjugate())                                                                                                           
 
stream_zs = np.zeros(Nx * Ny, complex)
for k in range(Nx * Ny):
    stream_zs[k] = stream(f_zs[k])
    
stream_xs = stream_zs.real
stream_ys = stream_zs.imag
    


#------------------------------------------------------------------------------
# Plot stream

title = "Stream Plot"
ax_labels = ["x", "y"]

# graph.plot_2D(stream_xs, f_xs, title, ["x", "$\phi_x(x,y)$"])
# graph.plot_2D(stream_xs, f_ys, title, ["y", "$\phi_x(x,y)$"])
# graph.plot_2D(stream_ys, f_xs, title, ["x", "$\phi_y(x,y)$"])
# graph.plot_2D(stream_ys, f_ys, title, ["y", "$\phi_y(x,y)$"])

zip_fxy = np.array(sorted(zip(f_xs, f_ys, stream_xs, stream_ys))).T  #sort in x
f_xs = zip_fxy[0]
f_ys = zip_fxy[1]
stream_xs = zip_fxy[2]
stream_ys = zip_fxy[3]

graphics.plot_stream(stream_xs, stream_ys, f_xs, f_ys, title, ax_labels)