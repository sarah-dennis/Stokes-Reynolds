# -*- coding: utf-8 -*-
"""
Created on Tue Oct  1 11:56:12 2024

link to sd_run.sh

@author: sarah
"""

import stokes_control as control
import stokes_examples as examples
#------------------------------------------------------------------------------
# Instructions: 
#------------------------------------------------------------------------------
 
# I. Uncomment and select an example. 
#    class structure: examples < stokes_examples < stokes_heights < domain
#           examples sets fluid parameters Re, Q, U, 
#           stokes_examples sets geometric parameters H, L, h, delta, etc. 
#           stokes_heights sets grid size N=1/dx and discretizes the height function
#        /> example = examples.BFS_H2L4_Re0_Q2_U0

# II. Look in ./examples/[example_name] for precomputed solutions to each example at various N

# III. Initialize the solver with the selected example
#        /> solver = control.Stokes_Solver(example)  

# IV. Choose a method from the solver, provide a grid size N to run or load 
#       the solver in stokes_control has several solution, graphing and convergence methods,
#     For basic use... 
#        /> solver.new_run(N) 
#        /> solver.load_run(N) 
#        /> solver.load_plot(N)

#------------------------------------------------------------------------------
# BFS Re=0
#------------------------------------------------------------------------------
# example = examples.BFS_H2p75L4_Re0_Q2_U0
# example = examples.BFS_H2p5L4_Re0_Q2_U0
# example = examples.BFS_H2p25L4_Re0_Q2_U0
# example = examples.BFS_H2L4_Re0_Q2_U0 # ***convergence N0=20,dN=2,Nmax=320****
# example = examples.BFS_H1p25L4_Re0_Q2_U0
# example = examples.BFS_H1p125L4_Re0_Q2_U0

example = examples.BFS_H2L4_Re0_Q1_U0
#------------------------------------------------------------------------------
# BFS Re=0.25
#------------------------------------------------------------------------------
# example = examples.BFS_H2p75L4_Re0p25_Q2_U0
# example = examples.BFS_H2p5L4_Re0p25_Q2_U0
# example = examples.BFS_H2p25L4_Re0p25_Q2_U0
# example = examples.BFS_H2L4_Re0p25_Q2_U0
# example = examples.BFS_H1p5L4_Re0p25_Q2_U0
# example = examples.BFS_H1p25L4_Re0p25_Q2_U0
# example = examples.BFS_H1p125L4_Re0p25_Q2_U0
#------------------------------------------------------------------------------
# BFS Re=0.5
#------------------------------------------------------------------------------
# example = examples.BFS_H2p75L4_Re0p5_Q2_U0
# example = examples.BFS_H2p5L4_Re0p5_Q2_U0
# example = examples.BFS_H2p25L4_Re0p5_Q2_U0
# example = examples.BFS_H2L4_Re0p5_Q2_U0
# example = examples.BFS_H1p5L4_Re0p5_Q2_U0
# example = examples.BFS_H1p25L4_Re0p5_Q2_U0
# example = examples.BFS_H1p125L4_Re0p5_Q2_U0
#------------------------------------------------------------------------------
# BFS Re=1
#------------------------------------------------------------------------------
# example = examples.BFS_H2p75L4_Re1_Q2_U0
# example = examples.BFS_H2p5L4_Re1_Q2_U0
# example = examples.BFS_H2p25L4_Re1_Q2_U0
# example = examples.BFS_H2L4_Re1_Q2_U0
# example = examples.BFS_H1p5L4_Re1_Q2_U0
# example = examples.BFS_H1p25L4_Re1_Q2_U0
# example = examples.BFS_H1p125L4_Re1_Q2_U0
#------------------------------------------------------------------------------
# smoothed-step H=2 Re=0
#------------------------------------------------------------------------------
# example = examples.dBFS_H2L4_d1p25_Re0_Q2_U0 
# example = examples.dBFS_H2L4_d1_Re0_Q2_U0
# example = examples.dBFS_H2L4_d0p75_Re0_Q2_U0
# example = examples.dBFS_H2L4_d0p5_Re0_Q2_U0 
# example = examples.dBFS_H2L4_d0p25_Re0_Q2_U0
# example = examples.dBFS_H2L4_d0p125_Re0_Q2_U0
# example = examples.dBFS_H2L4_d0_Re0_Q2_U0
#------------------------------------------------------------------------------
## smoothed-step H=1.5 Re=0
#------------------------------------------------------------------------------
# example = examples.dBFS_H1p5L4_d0p75_Re0_Q2_U0
# example = examples.dBFS_H1p5L4_d0p5_Re0_Q2_U0
# example = examples.dBFS_H1p5L4_d0p25_Re0_Q2_U0 # ***convergence N0=20,dN=2,Nmax=320****
# example = examples.dBFS_H1p5L4_d0p125_Re0_Q2_U0
# example = examples.dBFS_H1p5L4_d0_Re0_Q2_U0
#------------------------------------------------------------------------------
## smoothed-step H=1.25 Re=0
#------------------------------------------------------------------------------
# example = examples.dBFS_H1p25L4_d0_Re0_Q2_U0    
# example = examples.dBFS_H1p25L4_d0p05_Re0_Q2_U0                                      
# example = examples.dBFS_H1p25L4_d0p125_Re0_Q2_U0 
# example = examples.dBFS_H1p25L4_d0p25_Re0_Q2_U0
#------------------------------------------------------------------------------
## wedged corner BFS Re=0
#------------------------------------------------------------------------------
# example = examples.BFS_H2p75L4_noEddy_Re0_Q2_U0
# example = examples.BFS_H2p5L4_noEddy_Re0_Q2_U0
# example = examples.BFS_H2p25L4_noEddy_Re0_Q2_U0
# example = examples.BFS_H2L4_noEddy_Re0_Q2_U0
# example = examples.BFS_H1p5L4_noEddy_Re0_Q2_U0
# example = examples.BFS_H1p25L4_noEddy_Re0_Q2_U0
#------------------------------------------------------------------------------
## wedged corner variationds Re=0
#------------------------------------------------------------------------------
# example = examples.BFS_H2L4_cornerTriA_Re0_Q2_U0
# example = examples.BFS_H2L4_cornerTriC_Re0_Q2_U0
# example = examples.BFS_H2L4_cornerTriB_Re0_Q2_U0
# example = examples.BFS_H2L4_cornerTriD_Re0_Q2_U0
#------------------------------------------------------------------------------

# example = examples.TriCavity
# example = examples.bump_Re0
example=examples.logistic_Re0

#------------------------------------------------------------------------------

solver = control.Stokes_Solver(example, max_iters=5000)                

N=160

       
zoom_on=False               
# 
# solver.new_run(N) 

# solver.load_scale(20,40) 
# 
solver.load_run(N)

# 
# solver.load_run_new_many(N, 2, 3)  

solver.load_plot(N, zoom=zoom_on)

# ------------------------------------------------------------------------------
# solver.compare(20,[40,80,160],320)







