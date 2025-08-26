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

H=2
h=1
L=4
l=1

delta = 32 #slope: -delta*(H-h)/4


# args = [h, H, l, L]
# example = examples.BFS

# args = [H, L, delta]-
# example = examples.BFS_pwl

args = [H, h, L, delta]
example = examples.Logistic

# args = [H,h, L, delta]
# example = examples.Sinusoid

# args = [H, L]
# example = examples.TriCavity

# args =  [h, h*2, h, l, 1.25*l, 0.75*l, l]
# example = examples.TriSlider



# ------------------------------------------------------------------------------


solver = control.Stokes_Solver(example, args, U, Q, Re, max_iters=50000)                

N=160
       
zoom_on=  False               

# solver.new_run(N) 


# solver.load_scale(20,40) 
# solver.load_run(N)

# solver.load_run_many(20, 2, 2)

# solver.new_run_many(N, 2, 4)  
# solver.load_run_new_many(N, 2, 3)

solver.load_plot(N, zoom=zoom_on)

# ------------------------------------------------------------------------------
# solver.compare(20,[40,80,160],320)







