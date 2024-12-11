# -*- coding: utf-8 -*-
"""
Created on Thu Oct 24 13:23:27 2024

@author: sarah
"""

import reyn_control as control
import reyn_examples as examples

# Example = examples.BFS

Example = examples.BFS_deltaSmooth
# Example = examples.TriSlider
# Example = examples.BFS_noEddy
# Exyample = examples.HexSlider
# Example = examples.variableSlider
U= 0
dP = -61.75
N = 200


solver = control.Reynolds_Solver(Example, U, dP)
solver.solve_and_plot(N)

# control.convg_pwl_fd(Example, U, dP, N0=20, dN=2, many=6)
# 
