# -*- coding: utf-8 -*-
"""
Created on Tue Oct  1 11:56:12 2024

link to sd_run.sh

@author: sarah
"""

import stokes_control as control
import stokes_examples as examples

U=0
Q=1
Re=0

h_in = 2
h_out = 1
l_in = 8
l_out=8
args = [h_in, h_out, l_in, l_out]
Example = examples.BFS

# args = [H, L, delta]
# Example = examples.BFS_pwl

# H=2
# h=1
# L=16
# delta = 32 #slope: -delta*(H-h)/4
# args = [H, h, L, delta]
# Example = examples.Logistic

# args = [H,h, L, delta]
# Example = examples.Sinusoid

# args = [H, L]
# Example = examples.TriCavity

# l=7
# h=1
# H = 2
# args =  [h, H, h, l, 1.25, 0.75, l]
# Example = examples.TriSlider



# ------------------------------------------------------------------------------


solver = control.Stokes_Solver(Example, args, U, Q, Re, max_iters=500000)                

# N=80

zoom_on= not True           

# solver.new_run(20) 
# 

# solver.load_scale(80,160) 
solver.load_run(20)

# solver.load_run_many(20, 2, 4)

# solver.new_run_many(N, 2, 4)  
# solver.load_run_new_many(N, 2, 3)
# 
solver.load_plot(20, zoom=zoom_on)

# ------------------------------------------------------------------------------
# solver.compare(20,[40,80,160],320)







