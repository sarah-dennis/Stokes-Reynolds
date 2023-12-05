#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 30 11:05:44 2023

@author: sarahdennis
"""
import numpy as np
import graphics as graph

#------------------------------------------------------------------------------    
# t : Upper half plane
#             .
#             |
#             |
#             |
# -infty ---- 0 ---- 1 --- a ---- +infty
#  <-- L --- N=M --- P --- R ---- S -->

# t = x + iy

t_Nx = 20
t_Ny = 20

t_xMin = -10
t_xMax = 10
t_yMin = 0
t_yMax = 10

t_dx = (t_xMax - t_xMin)/(t_Nx-1)
t_dy = (t_yMax - t_yMin)/(t_Ny-1)

t_xs = np.zeros(t_Nx)
t_ys = np.zeros(t_Ny)

for i in range (t_Nx):
    t_xs[i] = t_xMin + i*t_dx
for j in range (t_Ny):
    t_ys[j] = t_yMin + j*t_dy

# t_plane = np.zeros((t_Ny, t_Nx), dtype=complex)
# for i in range(t_Nx):
#     for j in range(t_Ny):
#         t_plane[j][i] = complex(t_xs[i], t_ys[j])

#------------------------------------------------------------------------------
# Pentagon
#------------------------------------------------------------------------------
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

h1 = 1 #dist N -> M
h2 = 2 #dist S -> L

# Step height : 
#    h2 > h1, a = (H/h)**2
     
#------------------------------------------------------------------------------   
# # Schwarz-Christoffel Mapping z := f(t)   (z-plane)
#------------------------------------------------------------------------------       
a = (h2/h1)**2
b1 = h2/np.pi 
b2 = h1/np.pi
c = complex(0, -h2+h1)

def f(t):
    
    arg1 = (2*t-a-1)/(a-1)
    
    arg2 = ((a+1)*t - 2*a)/((a-1)*t)
    
    z = b1 * np.arccosh(arg1) - b2 * np.arccosh(arg2) + c
    
    return z

#---------------------------------------------------------------------------- 
z_Nx = 20
z_Ny = 20

z_xs = np.zeros(z_Nx)
z_ys = np.zeros(z_Ny)

z_xMin = -l1
z_xMax = l2
z_yMin = h1-h2
z_yMax = h2

z_dx = np.abs(z_xMin - z_xMax)/(z_Nx-1)
for i in range(z_Nx):
    z_xs[i] = z_xMin + i*z_dx
   
z_dy = np.abs(z_yMin - z_yMax)/(z_Ny-1)
for j in range(z_Ny):
   z_ys[j] = z_yMin + j*z_dy
   
z_plane = np.zeros((z_Ny, z_Nx), dtype=complex)

for i in range(z_Nx):
    for j in range(z_Ny):
        t = complex(t_xs[i], t_ys[j])
        z = complex(z_xs[i], z_ys[j])
        ft = f(t)
        z_plane[j][i]= ft
        print("z: (%.2f, %.2f), f(t): (%.2f, %.2f)" % (z.real, z.imag, ft.real, ft.imag))
# print(z_plane)
    
#------------------------------------------------------------------------------
# Stream Function Mapping   s := phi(z)
#------------------------------------------------------------------------------

d = complex(0, -U*h2/(2*np.pi))
def stream(z):
    return  d*np.log(z/z.conjugate())
   
s_x = np.zeros((z_Ny, z_Nx))
s_y = np.zeros((z_Ny, z_Nx))

# mask with height function...

for i in range(z_Nx):
    for j in range (z_Ny):
        if z_xs[i] >= 0 or z_ys[j] >= 0:
            s = stream(complex(z_xs[i], z_ys[j]))
            s_x[j][i] = s.real
            s_y[j][i] = s.imag


#------------------------------------------------------------------------------
# Plot stream

title = "Stream Plot"
ax_labels = ["x", "y"]

graph.plot_stream(s_x, s_y, z_xs, z_ys, title, ax_labels)

graph.plot_quivers(s_x, s_y, z_xs, z_ys, title, ax_labels)
    

    
    
    
    
    
    
    
    
    
    
    
    
    