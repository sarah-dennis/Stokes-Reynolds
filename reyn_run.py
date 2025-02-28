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
write_on=True
#------------------------------------------------------------------------------
## space parmaters |--> args = [ ... ]
#------------------------------------------------------------------------------
l = 0          # inlet/outlet length        
h = .125        # minimum height       
H = 4           # maximum height
xr = 0.35       # BFS reattachment point x=l+xr
yr = 0.41       # BFS detachment point y=yr
delta = 0.25    # smoothed step slope reciprocal 0 <= delta < L/2
r=1             # radius 0 <= r <= 1
k= 2* 3.14      # period sin(kx)

## args are specific to each example and domain is initialized at solve time

#------------------------------------------------------------------------------
## Piecewise-linear examples 
##       (analytic or finite difference solution)
#------------------------------------------------------------------------------
# Example = examples.BFS
# args = [H,l]

# Example = examples.BFS_noEddy
# args = [H,xr,yr] 

# Example = examples.BFS_deltaSmooth
# args = [H,delta]

# Example = examples.TriSlider
# args = None

# Example = examples.HexSlider
# args = None

# Example = examples.variableSlider  
# args = None

Example = examples.TriCavity
args = [H,l]

#------------------------------------------------------------------------------
## Smooth examples  
##      (finite difference solution only)
#------------------------------------------------------------------------------
# Example = examples.Bump # 
# args=[r,k]

Example = examples.CircCavity
args=[r,h,l]

#------------------------------------------------------------------------------
## surface velocity
U=0.5             # u(x,0) = U; u(x,h) = 0
                  # v(x,0) = 0; v(x,h) = 0

## pressure drop                 
dP = 0            # p(0,y) = -dP; p(xf,y) = 0


solver = control.Reynolds_Solver(Example, U, dP, args)
#------------------------------------------------------------------------------

## grid size
N = 160        # N = 1/dx = 1/dy

## solution methods (plots and returns pressure, velocity )
# solver.fd_solve(N, plot=plots_on, zoom=zoom_on)
# solver.pwl_solve(N, plot=plots_on, zoom=zoom_on)
solver.fd_adj_solve(N, write=write_on, plot=plots_on, zoom=zoom_on)

#------------------------------------------------------------------------------
# solver.load_plot(N, zoom=zoom_on)

#------------------------------------------------------------------------------
# solver.convg_pwl_fd([20,40,80,160])
# solver.convg_adj_fd(20, [40, 80, 160, 320,640,1280],2560)