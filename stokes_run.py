# -*- coding: utf-8 -*-
"""
Created on Tue Oct  1 11:56:12 2024

link to sd_run.sh

@author: sarah
"""

import stokes_control as control
import stokes_examples as examples
# 
# example = examples.BFS_H2L4_Re0_Q2_U0

# example = examples.BFS_H2L4_delta1p5 
# example = examples.BFS_H2L4_delta1
# example = examples.BFS_H2L4_delta0p75
# example = examples.BFS_H2L4_delta0p5 
# example = examples.BFS_H2L4_delta0p25
#...
# example = examples.BFS_H2L4_delta0


# example = examples.BFS_H1p5L4_Re0_Q2_U0


# example = examples.BFS_H1p5L4_delta1
# example = examples.BFS_H1p5L4_delta0p75
# example = examples.BFS_H1p5L4_delta0p5 
# example = examples.BFS_H1p5L4_delta0p25
example = examples.BFS_H1p5L4_delta0p125
# example = examples.BFS_H1p5L4_delta0


# example = examples.BFS_H1p25L4_Re0_Q2_U0


# example = examples.BFS_H1p125L4_Re0_Q2_U0

N = 200
#------------------------------------------------------------------------------
solver = control.Stokes_Solver(example) 
#------------------------------------------------------------------------------
# solver methods: new_run, load_run, new_run_many, load_run_new_many, load_run_many
#------------------------------------------------------------------------------
# solver.new_run(200, solver.max_iters)
solver.load_run_new_many(N, 2, 1)
#------------------------------------------------------------------------------
# solver.load_plot(N, zoom=True)
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
# BFS EXAMPLES  
#------------------------------------------------------------------------------
# BFS H=2, L=4, Re=0, Q=2, U=0 --> done
# solver.compare(20,[40,80,160,320],640)


# BFS H=1.5, L= 4 Re=0, Q=2, U=0 --> done
# solver.compare(20,[40,80,160,320],640)

# TODO: BFS H=1.25, L= 4 Re=0, Q=2, U=0 --> running N=640
# solver.compare(20,[40,80,160],320)


# TODO: BFS H=1.125, L= 4 Re=0, Q=2, U=0 --> N=200 running cluster

#------------------------------------------------------------------------------
# BFS EDDIE REMOVAL EXAMPLES
#------------------------------------------------------------------------------
# TODO: determine reattachment length from BFS examples
#------------------------------------------------------------------------------
# Delta-BFS Examples
#------------------------------------------------------------------------------
# H=2, 1.5, 1.25, 1.125 delta=1.5, 1, 0.75, 0.5, 0.25, 0.125, 0

# dBFS d=0.75, H=2,L=4,Re=0,Q=2,U=0 --> done 
# solver.compare(20,[40,80,160,320],640)

# TODO: dBFS d=0.125, H=2,L=4,Re=0,Q=2,U=0 --> N=200 running cluster