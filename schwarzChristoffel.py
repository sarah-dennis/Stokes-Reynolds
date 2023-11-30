#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 30 11:05:44 2023

@author: sarahdennis
"""
import numpy as np
import graphics as graph


#------------------------------------------------------------------------------
# Pentagon

# -  M ----------- L -
# h1 |             | |
# -  N ---- P      | h2
#           |      | |
#           R ---- S -     
#    |--l1--|--l2--|  

# flow:  <--U-- 
U=0.5

l1 = 2 # dist N -> P
l2 = 5 # dist R -> S

h1 = 2 #dist N -> M
h2 = 3 #dist S -> L

# Step height : 
#    h2 > h1, a = (H/h)**2

#------------------------------------------------------------------------------
# Conformally Transformed Plane     (t-plane)

#              |      *
#           *  |              *
#       *      |                   *
#    *         |                       *
#              |    
# L --------- M~N ---- P ------- R ----- S
# t0 -(l1+l2)  0   +l1    +h2-h1    +l2 trad

Nx = 50
Ny = 50

t_x0 = -(l1+l2)
t_y_0 = 0

t_xlen = 2*(l1+l2) + h2-h1
t_ylen = t_xlen/2

t_dx = t_xlen / (Nx-1)
t_dy = t_ylen/ (Ny-1)

t_plane = np.zeros((Ny, Nx), dtype=complex)

for i in range(Nx):
    x = t_x0 + i*t_dx
    
    for j in range(Ny):
        y = j*t_dy
        
        t_plane[j][i] = complex(x,y)

        
#------------------------------------------------------------------------------    
# Schwarz-Christoffel Mapping z := f(t)   (z-plane)
        
a = (h2/h1)**2
b1 = h2/np.pi 
b2 = h1/np.pi
c = complex(0,-np.pi*(1-h1/h2))

def f(t):
    arg1 = (2*t - (a+1))/(a-1)
    
    arg2 = ((a+1)*t - 2*a)/((a-1)*t)
    
    z = b1 * np.arccosh(arg1) - b2 * np.arccosh(arg2) + c
    
    return z
    
    
z_plane = np.zeros((Ny, Nx), dtype=complex)

for i in range(Nx):
    for j in range (Ny):
        z_plane[j][i]=f(t_plane[j][i])
    
#------------------------------------------------------------------------------
# Stream Function Mapping   s := phi(z)


# --- for plotting make C grid ----
xs = np.zeros(Nx)
ys = np.zeros(Ny)

z_x0 = -l1
z_y0 = h1-h2

z_dx = (l1 + l2)/(Nx-1)
for i in range(Nx):
    xs[i] = z_x0 + i*z_dx
   
z_dy = h2/(Ny-1)
for j in range(Ny):
   ys[j] = z_y0 + j*z_dy
    
#-----------------------------------
# compute stream s = phi(z) 


d = complex(0, U*h2/(2*np.pi))
def stream(z):
       return d * np.log(z/z.conjugate())
   
s_x = np.zeros((Ny, Nx))
s_y = np.zeros((Ny, Nx))

# mask with height function...

for i in range(Nx):
    for j in range (Ny):
        if xs[i] >= 0 or ys[j] >= 0:
        
            s = stream(z_plane[j][i])
            s_x[j][i] = s.real
            s_y[j][i] = s.imag


#------------------------------------------------------------------------------
# Plot stream

title = "Stream Plot"
ax_labels = ["x", "y"]


# graph.plot_stream(s_x, s_y, xs, ys, title, ax_labels)

graph.plot_quivers(s_x, s_y, xs, ys, title, ax_labels)
    

    
    
    
    
    
    
    
    
    
    
    
    
    