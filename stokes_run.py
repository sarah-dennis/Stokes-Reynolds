# -*- coding: utf-8 -*-
"""
Created on Tue Oct  1 11:56:12 2024

link to sd_run.sh

@author: sarah
"""

import stokes_control as control
import stokes_examples as examples

#------------------------------------------------------------------------------
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
example = examples.BFS_H2L4_delta0p125
# example = examples.BFS_H2L4_delta0
#------------------------------------------------------------------------------
# example = examples.BFS_H1p5L4_delta1
# example = examples.BFS_H1p5L4_delta0p75
# example = examples.BFS_H1p5L4_delta0p5 
# example = examples.BFS_H1p5L4_delta0p25
# example = examples.BFS_H1p5L4_delta0p125 
# example = examples.BFS_H1p5L4_delta0
#------------------------------------------------------------------------------
# example = examples.BFS_H1p25L4_delta0                                             
# example = examples.BFS_H1p25L4_delta0p125 
#------------------------------------------------------------------------------
# example = examples.BFS_H2p75L4_noEddy_Re0_Q2_U0
# example = examples.BFS_H2p5L4_noEddy_Re0_Q2_U0
# example = examples.BFS_H2p25L4_noEddy_Re0_Q2_U0
# example = examples.BFS_H2L4_noEddy_Re0_Q2_U0
# example = examples.BFS_H1p5L4_noEddy_Re0_Q2_U0
# example = examples.BFS_H1p25L4_noEddy_Re0_Q2_U0

# example = examples.BFS_H2L4_cornerTriA_Re0_Q2_U0
# example = examples.BFS_H2L4_cornerTriB_Re0_Q2_U0
# example = examples.BFS_H2L4_cornerTriC_Re0_Q2_U0
# example = examples.BFS_H2L4_cornerTriD_Re0_Q2_U0
# example = examples.BFS_H2L4_cornerTriE_Re0_Q2_U0



# example = examples.basic
# example = examples.BFS_biswas_Re0

N=200
#------------------------------------------------------------------------------
solver = control.Stokes_Solver(example)                                      

# solver.new_run(N, solver.max_iters) 
# solver.load_scale(100,N)

solver.load_run(N,1)

# solver.load_run(N,solver.max_iters)

solver.load_plot(N, zoom=False)
# solver.load_plot(N, zoom=True)
# ------------------------------------------------------------------------------
