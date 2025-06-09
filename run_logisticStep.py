# -*- coding: utf-8 -*-
"""
Created on Thu May 29 17:16:18 2025

@author: sarah
"""

import stokes_control as control
import stokes_examples as examples

#--> DONE
# example = examples.logisticBFS_l1_Re0_Q1p19_U0
# example = examples.logisticBFS_l0p5_Re0_Q0p99_U0
# example = examples.logisticBFS_l0p25_Re0_Q0p85_U0
# example = examples.logisticBFS_l0p125_Re0_Q0p79_U0

#------------------------------------------------------------------------------
# --> new examples, start at N=20

#------------------------------------------------------------------------------
example = examples.logisticBFS_l0p0625_Re0_Q0p77_U0

solver = control.Stokes_Solver(example, max_iters=50000)          
N=160
solver.load_run(N) 

#------------------------------------------------------------------------------
example = examples.logisticBFS_l0p03125_Re0_Q0p75_U0

solver = control.Stokes_Solver(example, max_iters=50000)          
N_0=20
dN=2 
many=4 #(20, 40, 80, 160, 320)
solver.new_run_many(N_0, dN, many) 
