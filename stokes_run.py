# -*- coding: utf-8 -*-
"""
Created on Tue Oct  1 11:56:12 2024

link to sd_run.sh

@author: sarah
"""

import stokes_control as control
import stokes_examples as examples


zoom_on= not True    

U=0
Q=1
Re=0


#------------------------------------------------------------------------------
# h_in = 2.75
# h_out = 1
# l_in = 8
# l_out=8
# args = [h_in, h_out, l_in, l_out]
# Example = examples.BFS


#------------------------------------------------------------------------------
# h_in = 2
# h_out = 1
# l_in = 8
# l_out=8

# xr = 0.35
# yr = 0.4
# xr_0p75 = 0.2625
# yr_0p75 = 0.3
# xr_0p5 = 0.175
# yr_0p5 = 0.2

# args = [h_in, h_out, l_in, l_out, xr, yr]
# Example = examples.BFS_wedge


# args = [h_in, h_out, l_in, l_out, xr_0p5, yr_0p5]
# Example = examples.BFS_wedge

# args = [h_in, h_out, l_in, l_out, xr_0p75, yr_0p75]
# Example = examples.BFS_wedge

H=2
h=1
L=16
delta=1/4
args = [H, h, L, delta]
Example = examples.BFS_pwl

#------------------------------------------------------------------------------
H=2
h=1
L=16
delta=1/4
args = [H, h, L, delta]
Example = examples.BFS_pwl


#------------------------------------------------------------------------------
# H=2
# h=1
# L=4
# delta = 2 #slope: -delta*(H-h)/4
# args = [H, h, L, delta]
# Example = examples.Logistic


#------------------------------------------------------------------------------
# H = 2
# delta = 1/4
# k = 2 #int 
# l = 1
# L = 8
# args = [H, delta, k,  l, L]
# Example = examples.Sinusoid


#------------------------------------------------------------------------------
# H = 2
# L = 2 
# args = [H, L]# tri slope = 2H/L
# Example = examples.TriCavity


#------------------------------------------------------------------------------
# l=7
# h=1
# H = 2
# args =  [h, H, h, l, 1.25, 0.75, l]
# Example = examples.TriSlider



#------------------------------------------------------------------------------
solver = control.Stokes_Solver(Example, args, U, Q, Re, max_iters=500000)                

N=40

       

# solver.new_run(N) 
# solver.load_run(N)

# solver.load_scale(N,2*N) 

# solver.load_copy(N, new_Example, new_args)

# solver.load_run_many(N, 2, 4)

solver.new_run_many(N, 2, 3)  
# solver.load_run_new_many(N, 2, 2)

solver.load_plot(N, zoom=zoom_on)

# ------------------------------------------------------------------------------
# solver.compare(20,[40,80,160],320)







