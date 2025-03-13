# -*- coding: utf-8 -*-
"""
Created on Thu Oct 24 13:23:27 2024

@author: sarah
"""

import reyn_control as control
import reyn_examples as examples

#------------------------------------------------------------------------------
## Instructions: 
#------------------------------------------------------------------------------

# I. Set up the example. 
#       for examples, see reyn_examples, and reyn_heights superclass
#       set args (geometric parameters) H, L, h, r, etc.
#       
#       /> Example = examples.BFS
#       /> args = [H=2,l=4] 

# II. Initialize the solver 
#       set the boundary conditions (U, dP)

#        /> U, dP = [100, 1 , 0]
#        /> solver = control.Reynolds_Solver(Example, U, dP, args)

# III. Choose a method for the solver, provide a grid size N to run or load 
#       the solver in reyn_control has several solution, graphing and convergence methods,

#        /> solver.fd_solve(N, plot=plots_on)    # Finite difference solve

#        /> solver.pwl_solve(N, plot=plots_on)   # Analytic solve for piecewise-linear heights

#        /> solver.fd_adj_solve(N, plot=plots_on) #Finite difference solve for *adjusted Reynolds equation*
#------------------------------------------------------------------------------
plots_on=True
zoom_on=False 
write_on=False
#------------------------------------------------------------------------------
## space parmaters |--> args = [ ... ]
#------------------------------------------------------------------------------
# l = 1     # inlet/outlet length        
# h = .125        # minimum height       
# H = 1/2           # maximum height
# xr = 0.35       # BFS reattachment point x=l+xr
# yr = 0.41       # BFS detachment point y=yr
# delta = 0.5    # smoothed step slope reciprocal 0 <= delta < L/2
# r = 0.             # radius 0 <= r < 1/2
# lam = 0.2
# k= 2      # period sin(kx)
# L= 4
## args are specific to each example and domain is initialized at solve time

#------------------------------------------------------------------------------
## Piecewise-linear examples 
##       (analytic or finite difference solution)
#------------------------------------------------------------------------------
# Example = examples.BFS
# args = [H,l, L]

# Example = examples.BFS_noEddy
# args = [H,xr,yr] 

# Example = examples.BFS_deltaSmooth
# args = [H,delta]

# Example = examples.TriSlider
# args = None

# Example = examples.HexSlider
# args = None

Example = examples.TriCavity
H=4
L=1
args = [H,L]

#------------------------------------------------------------------------------
## Smooth examples  
##      (finite difference solution only)
#------------------------------------------------------------------------------
# Example = examples.Sinusoid
# args = [r, k, L]

# Example = examples.LambdaBump # 
# lam=0.2
# H=1
# l=1
# args=[lam, H, l]

# Example = examples.Cylinder
# r=1
# h0 = 1
# l=2
# #d = h0+r
# args= [r,h0,l]

U=1
dP=0

#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
# solution methods (plots  and returns pressure, velocity )

N=100

solver = control.Reynolds_Solver(Example, U, dP, args)
# solver.fd_solve(N, plot=plots_on, zoom=zoom_on)
solver.pwl_solve(N, plot=plots_on, zoom=zoom_on)
solver.fd_adj_solve(N, write=write_on, plot=plots_on, zoom=zoom_on)
solver.fd_pert_solve(N, write=write_on, plot=plots_on, zoom=zoom_on)

#------------------------------------------------------------------------------
# solver.load_plot(N, zoom=zoom_on)

#------------------------------------------------------------------------------
# solver.convg_pwl_fd([20,40,80,160])

