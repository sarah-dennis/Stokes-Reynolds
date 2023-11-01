# -*- coding: utf-8 -*-
"""
Created on Tue Oct 31 15:07:45 2023

BFS hydraulic resistance
Stokes vs Reynolds comparison

@author: sarah
"""

import domain as dm

# x boundary
x0 = 0

xf = 5    
BC = "fixed" 

# surface velocity 
U = 1   #V_x(x,0) = U

# kinematic viscosity
visc = 1  

#friction coeff (Stokes)
lam = 48 

# Grid size
Nx = 500

domain = dm.Domain(x0, xf, visc, U, Nx, BC)

h1 = 1 
h2 = 1.9423*h1

#bulk velocity
Ub1 = 3 #...?
Ub2 = Ub1 * h1/h2



