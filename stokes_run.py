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
# example = examples.BFS_H2L4_Re0_Q2_U0 # ***convergence N0=20,dN=2,Nmax=320****
# example = examples.BFS_H1p25L4_Re0_Q2_U0
# example = examples.BFS_H1p125L4_Re0_Q2_U0
#------------------------------------------------------------------------------
# example = examples.BFS_H2p75L4_Re0p25_Q2_U0
# example = examples.BFS_H2p5L4_Re0p25_Q2_U0
# example = examples.BFS_H2p25L4_Re0p25_Q2_U0
# example = examples.BFS_H2L4_Re0p25_Q2_U0
# example = examples.BFS_H1p5L4_Re0p25_Q2_U0
# example = examples.BFS_H1p25L4_Re0p25_Q2_U0
# example = examples.BFS_H1p125L4_Re0p25_Q2_U0
#------------------------------------------------------------------------------
# example = examples.BFS_H2p75L4_Re0p5_Q2_U0
# example = examples.BFS_H2p5L4_Re0p5_Q2_U0
# example = examples.BFS_H2p25L4_Re0p5_Q2_U0
# example = examples.BFS_H2L4_Re0p5_Q2_U0
# example = examples.BFS_H1p5L4_Re0p5_Q2_U0
# example = examples.BFS_H1p25L4_Re0p5_Q2_U0
# example = examples.BFS_H1p125L4_Re0p5_Q2_U0
#------------------------------------------------------------------------------
# example = examples.BFS_H2p75L4_Re1_Q2_U0
# example = examples.BFS_H2p5L4_Re1_Q2_U0
# example = examples.BFS_H2p25L4_Re1_Q2_U0
# example = examples.BFS_H2L4_Re1_Q2_U0
# example = examples.BFS_H1p5L4_Re1_Q2_U0
# example = examples.BFS_H1p25L4_Re1_Q2_U0
# example = examples.BFS_H1p125L4_Re1_Q2_U0
#------------------------------------------------------------------------------
# example = examples.dBFS_H2L4_d1p25_Re0_Q2_U0 
# example = examples.dBFS_H2L4_d1_Re0_Q2_U0
# example = examples.dBFS_H2L4_d0p75_Re0_Q2_U0
# example = examples.dBFS_H2L4_d0p5_Re0_Q2_U0 
# example = examples.dBFS_H2L4_d0p25_Re0_Q2_U0
# example = examples.dBFS_H2L4_d0p125_Re0_Q2_U0
# example = examples.dBFS_H2L4_d0_Re0_Q2_U0
#------------------------------------------------------------------------------
# example = examples.dBFS_H1p5L4_d0p75_Re0_Q2_U0
# example = examples.dBFS_H1p5L4_d0p5_Re0_Q2_U0
example = examples.dBFS_H1p5L4_d0p25_Re0_Q2_U0 # ***convergence N0=20,dN=2,Nmax=320****
# example = examples.dBFS_H1p5L4_d0p125_Re0_Q2_U0
# example = examples.dBFS_H1p5L4_d0_Re0_Q2_U0
#------------------------------------------------------------------------------
# example = examples.dBFS_H1p25L4_d0_Re0_Q2_U0    
# example = examples.dBFS_H1p25L4_d0p05_Re0_Q2_U0                                      
# example = examples.dBFS_H1p25L4_d0p125_Re0_Q2_U0 
# example = examples.dBFS_H1p25L4_d0p25_Re0_Q2_U0
#------------------------------------------------------------------------------
# example = examples.BFS_H2p75L4_noEddy_Re0_Q2_U0
# example = examples.BFS_H2p5L4_noEddy_Re0_Q2_U0
# example = examples.BFS_H2p25L4_noEddy_Re0_Q2_U0
# example = examples.BFS_H2L4_noEddy_Re0_Q2_U0
# example = examples.BFS_H1p5L4_noEddy_Re0_Q2_U0
# example = examples.BFS_H1p25L4_noEddy_Re0_Q2_U0
# ------------------------------------------------------------------------------
# example = examples.BFS_H2L4_cornerTriA_Re0_Q2_U0
# example = examples.BFS_H2L4_cornerTriC_Re0_Q2_U0
# example = examples.BFS_H2L4_cornerTriB_Re0_Q2_U0

# example = examples.BFS_H2L4_cornerTriD_Re0_Q2_U0



# example = examples.basic
# example = examples.BFS_biswas_Re0
N=200
#------------------------------------------------------------------------------
solver = control.Stokes_Solver(example)                                      

# solver.new_run(N, solver.max_iters) 
# solver.load_scale(100,200) 

# solver.load_run(N,1)

# solver.load_run(N,solver.max_iters)

# ## solver.load_run_new_many(N, 2, 2)

# solver.load_plot(N, zoom=False)
# print(solver.get_attachments(N))
# solver.load_plot(N, zoom=True)
 
# ------------------------------------------------------------------------------
# solver.compare(20,[40,80,160],320, linthresh=1e-8)







