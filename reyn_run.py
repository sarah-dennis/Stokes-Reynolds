# -*- coding: utf-8 -*-
"""
Created on Thu Oct 24 13:23:27 2024

@author: sarah
"""

import reyn_control as control
import reyn_examples as examples

# example = examples.BFS
example = examples.BFS_deltaSmooth
# example = examples.TriSlider
# example = examples.BFS_noEddy
# example = examples.HexSlider

U = 2
dP = -1
N = 100
solver = control.Reynolds_Solver(example, U, dP)
solver.solve_plot(N)
# solver.cnvg(N,2,6)