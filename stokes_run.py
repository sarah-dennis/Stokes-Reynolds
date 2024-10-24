# -*- coding: utf-8 -*-
"""
Created on Tue Oct  1 11:56:12 2024

link to sd_run.sh

@author: sarah
"""

import stokes_control as control
import stokes_examples as examples

example = examples.BFS_H2L4_Re0_Q2_U0

# example = examples.BFS_delta1p5 #N20
# example = examples.BFS_deltay1p0 #N20
# example = examples.BFS_delta0p75 #N20
# example = examples.BFS_delta0p5 #N20
# example = examples.BFS_delta0p25 #N20

# example = examples.BFS_noStep_Re0_Q2_U0
# example = examples.BFS_noEddy_Re0_Q2_U0
# example = examples.BFS_twostep02_Re0_Q2_U0

# example = examples.HexSlider_Re0_Q2_U0
# example = examples.HexSlider_Re05_Q2_U0
# example = examples.HexSlider_Re1_Q2_U0

N = 320
#------------------------------------------------------------------------------
solver = control.Stokes_Solver(example)
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# solver.new_run(N, 200)

# solver.load_run(320,solver.max_iters)

# solver.load_scale(320, 640)

# solver.load_run(640,solver.max_iters)
# solver.load_run_many(N,2,2)
#------------------------------------------------------------------------------
solver.load_plot(N)
#------------------------------------------------------------------------------
# HexSlider Re=0.5, Q=2, U=0
# solver.compare(20,[40,80,160,320],640)

# BFS H2L4 Re=0, Q=2, U=0
solver.compare(10,[20,40,80],160)

#BFS H2L4 delta=0.75, Re=0,Q=2,U=0 --> order 2??
# solver.compare(20,[40,80,160],320)
#------------------------------------------------------------------------------