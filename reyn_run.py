# -*- coding: utf-8 -*-
"""
Created on Thu Oct 24 13:23:27 2024

@author: sarah
"""

import reyn_control as control
import reyn_examples as examples

U= 0
dP = -54
N = 200
linthresh=1e-5
H=2
l=2
xr = 0.35 
yr = 0.41

delta = 0.25

Example = examples.BFS
args = [H,l]

# Example = examples.BFS_noEddy
# args = [H,xr,yr] 

# Example = examples.BFS_deltaSmooth
# args = [H,delta]


# Example = examples.TriSlider
# Exyample = examples.HexSlider
# Example = examples.variableSlider

solver = control.Reynolds_Solver(Example, U, dP, args)
flux = solver.solve_and_plot(N)

# control.convg_pwl_fd(Example, U, dP, args, N0=20, dN=2, many=4,linthresh=linthresh)
# 
