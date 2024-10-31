# -*- coding: utf-8 -*-
"""
Created on Thu Oct 24 13:23:27 2024

@author: sarah
"""

import reyn_control as control
import reyn_examples as examples

Example = examples.BFS
# Example = examples.BFS_deltaSmooth
# example = examples.TriSlider
# example = examples.BFS_noEddy
# example = examples.HexSlider

U = 0
dP = -76.09
N = 100

solver = control.Reynolds_Solver(Example, U, dP)
solver.solve_and_plot(N)

# control.convg_pwl_fd(Example, U, dP, N0=100, dN=2, many=5)

