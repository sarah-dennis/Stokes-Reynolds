# -*- coding: utf-8 -*-
"""
Created on Thu Oct 24 13:23:27 2024

@author: sarah
"""

import reyn_control as control
import reyn_examples as examples
plots_on=True

linthresh=1e-5


N = 300

U= 1
dP = 0


H=2
l=1

xr = 0.35 
yr = 0.41

delta = 0.5

r=0.5
k=3.14

#------------------------------------------------------------------------------
Example = examples.BFS
args = [H,l]

# Example = examples.BFS_noEddy
# args = [H,xr,yr] 

Example = examples.BFS_deltaSmooth
args = [H,delta]

# Example = examples.Bump #  #TODO currently fd only - not pwl 
# args=[r,k]

# Example = examples.TriCavity
# args = [H]

# Example = examples.TriSlider
# args = None

# Example = examples.HexSlider
# args = None

# Example = examples.variableSlider  
# args = None
#------------------------------------------------------------------------------

solver = control.Reynolds_Solver(Example, U, dP, args)
# solver.fd_solve(N, plot=plots_on)
# solver.pwl_solve(N, plot=plots_on)
solver.fd_adj_solve(N, plot=plots_on)

# control.convg_pwl_fd(Example, U, dP, args, N0=20, dN=2, many=4,linthresh=linthresh)
