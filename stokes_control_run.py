# -*- coding: utf-8 -*-
"""
Created on Tue Oct  1 11:56:12 2024

link to sd_run.sh

@author: sarah
"""

import stokes_control as control
import stokes_examples as examples

# example = examples.BFS_H2L4_Re0_Q2_U0

example = examples.BFS_smoothd0125

# example = examples.BFS_noStep_Re0_Q2_U0
# example = examples.BFS_noEddy_Re0_Q2_U0
# example = examples.BFS_twostep02_Re0_Q2_U0

# example = examples.HexSlider_Re0_Q2_U0
# example = examples.HexSlider_Re05_Q2_U0
# example = examples.HexSlider_Re1_Q2_U0
#------------------------------------------------------------------------------
solver = control.Stokes_Solver(example)
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------

# solver.new_run(200,solver.max_iters)
solver.load_run(200,solver.max_iters)



#------------------------------------------------------------------------------
solver.load_plot(200)

#------------------------------------------------------------------------------
# HexSlider Re=0.5, Q=2, U=0
# solver.compare(20,[40,80,160],320) 

#------------------------------------------------------------------------------