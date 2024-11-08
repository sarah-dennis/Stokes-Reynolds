# -*- coding: utf-8 -*-
"""
Created on Tue Oct  1 11:56:12 2024

link to sd_run.sh

@author: sarah
"""

import stokes_control as control
import stokes_examples as examples
import graphics
import numpy as np

# example = examples.BFS_H2p75L4_Re0_Q2_U0
# example = examples.BFS_H2p5L4_Re0_Q2_U0
# example = examples.BFS_H2p25L4_Re0_Q2_U0
# example = examples.BFS_H2L4_Re0_Q2_U0
# example = examples.BFS_H1p5L4_Re0_Q2_U0
# example = examples.BFS_H1p25L4_Re0_Q2_U0
# example = examples.BFS_H1p125L4_Re0_Q2_U0
#------------------------------------------------------------------------------
# example = examples.BFS_H2L4_delta1p5 
# example = examples.BFS_H2L4_delta1p25 
# example = examples.BFS_H2L4_delta1
# example = examples.BFS_H2L4_delta0p75 
# example = examples.BFS_H2L4_delta0p5 
# example = examples.BFS_H2L4_delta0p25
# example = examples.BFS_H2L4_delta0p125
example = examples.BFS_H2L4_delta0
#------------------------------------------------------------------------------
# example = examples.BFS_H1p5L4_delta0p75
# example = examples.BFS_H1p5L4_delta0p5 
# example = examples.BFS_H1p5L4_delta0p25
# example = examples.BFS_H1p5L4_delta0p125 
# example = examples.BFS_H1p5L4_delta0
#------------------------------------------------------------------------------
# example = examples.BFS_H2p5L4_delta1p5
# example = examples.BFS_H2p5L4_delta0p75
# example = examples.BFS_H2p5L4_delta0p5
# example = examples.BFS_H2p5L4_delta0p375
#------------------------------------------------------------------------------

# example = examples.BFS_H2p5L4_noEddy_Re0_Q2_U0
#------------------------------------------------------------------------------
# ex1 = examples.BFS_H2p5L4_noEddy_Re0_Q2_U0
# ex2 = examples.BFS_H2p25L4_noEddy_Re0_Q2_U0
# ex3 = examples.BFS_H2L4_noEddy_Re0_Q2_U0
# ex4 = examples.BFS_H1p5L4_noEddy_Re0_Q2_U0
# ex5 = examples.BFS_H1p25L4_noEddy_Re0_Q2_U0
# examples = [ex1,ex2,ex3,ex4,ex5]
#------------------------------------------------------------------------------

N = 40
# for example in examples:

#     solver = control.Stokes_Solver(example) 
#     solver.new_run(N,solver.max_iters)
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
solver = control.Stokes_Solver(example)
solver.load_run(N,solver.max_iters)
# solver.load_scale(N,2*N) 
# solver.load_run(2*N,solver.max_iters)
solver.load_plot (N, zoom=True)
# ------------------------------------------------------------------------------
