# -*- coding: utf-8 -*-
"""
Created on Mon Sep 22 12:15:51 2025

@author: sarah
"""

import stokes_control as control
import stokes_examples as examples

U=0
Q=1
Re=0


h=1
l=7
peak = 1/16

args =  [h, h*peak, h, l, l*1.25, l*0.75, l]
example = examples.TriSlider



# ------------------------------------------------------------------------------


solver = control.Stokes_Solver(example, args, U, Q, Re, max_iters=50000)                

N=80
solver.load_run(N)
    # ------------------------------------------------------------------------------
  
peaks =  [1/8, 1/4, 1/2, 3/4, 5/4, 3/2, 7/4, 2]

for peak in peaks:
    args =  [h, h*peak, h, l, l*1.25, l*0.75, l]
    example = examples.TriSlider
    
    
        
    
    solver = control.Stokes_Solver(example, args, U, Q, Re, max_iters=50000)                
    
    N=20
          
    solver.new_run_many(N, 2, 3)  
