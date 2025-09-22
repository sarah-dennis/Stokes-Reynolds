# -*- coding: utf-8 -*-
"""
Created on Mon Sep 22 11:50:47 2025

@author: sarah
"""


import stokes_control as control
import stokes_examples as examples

U=0
Q=1
Re=0

H=2
h=1
L=16
l=1
# ------------------------------------------------------------------------------
delta = 3

args = [H, h, L, delta]
example = examples.Logistic





solver = control.Stokes_Solver(example, args, U, Q, Re, max_iters=50000)                

N=80
  

solver.load_run(N) 
 # ------------------------------------------------------------------------------
deltas = [4,6,8,16,32] #slope: -delta*(H-h)/4

for delta in deltas:
    args = [H, h, L, delta]
    example = examples.Logistic



    solver = control.Stokes_Solver(example, args, U, Q, Re, max_iters=50000)                

    N=20
      

    solver.new_run_many(N, 2, 3)  




