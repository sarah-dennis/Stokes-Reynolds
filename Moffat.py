# -*- coding: utf-8 -*-
"""
Created on Tue Nov 14 15:43:27 2023

@author: sarah
"""
import domain 
import heights 

x0 = 0
xf = 1
BC = "fixed"

U = 1

visc = 1

Nx = 100

domain = domain.Domain(x0, xf, visc, U, Nx, BC)
 
p0 = 0
pN = 1

height = heights.