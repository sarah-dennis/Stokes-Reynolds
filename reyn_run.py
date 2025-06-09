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
plots_on = True
zoom_on = False 
write_on = False
#------------------------------------------------------------------------------
## Piecewise-linear examples 
##       (analytic or finite difference solution)
#------------------------------------------------------------------------------
Example = examples.BFS
H=2
l=1
L=4
args = [H, l, L]

# Example = examples.BFS_noEddy
# H = 2
# xr = .6
# yr = .4
# args = [H,xr,yr] 

# Example = examples.BFS_deltaSmooth
# H = 2
# delta = 0.5
# args = [H,delta]

# Example = examples.TriSlider
# args = None

# Example = examples.HexSlider
# args = None

# Example = examples.TriCavity
# H=4
# L=1
# h=1/10
# l=1
# args = [H,L,h,l]

#------------------------------------------------------------------------------
## Smooth examples  
##      (finite difference solution only)
#------------------------------------------------------------------------------
# Example = examples.Sinusoid
# r=0.8
# k=3.14
# L=2
# args = [r, k, L]

# Example = examples.LambdaBump 


# lam=-1/4
# H=1/2
# l=1
# h0 = 1/2
# args=[lam, H, l, h0] 



# Example = examples.Cylinder
# drdx = 0
# r=1
# h0 = 1/2
# l=0.05
# args= [ r, h0,l, drdx]

#d/a = 1.5
#dP = 0
#U0 = 1
#h0 = 1

# k = 1/8 # slope -1/k at x=0
# H = 1  # logistic height scale

# h = 1  # min height 
# l = 2  # half length

# args = [k,H,h, l]

#------------------------------------------------------------------------------
# boundary conditions
#------------------------------------------------------------------------------

# U: velocity {u(x,y0)=U, u(x,h(x))=0}  {v(x,y0)=0, v(x,h(x))=0} 
U =-1
# dP: 1D pressure {p(x0,y)=, u(x,h(x))=0} 

dP =0

#------------------------------------------------------------------------------
# solution methods (plots  and returns pressure, velocity )

N = 100

solver = control.Reynolds_Solver(Example, U, dP, args)
# solver.fd_solve(N, plot=plots_on, zoom=zoom_on, uv=True)
solver.pwl_solve(N, plot=plots_on, zoom=zoom_on)
solver.fd_adj_solve(N, write=write_on, plot=plots_on, zoom=zoom_on, uv=True)
solver.fd_pert_solve(N, order=4, write=write_on, plot=plots_on, zoom=zoom_on, uv=True)

#------------------------------------------------------------------------------
# solver.load_plot(N, zoom=zoom_on)

#------------------------------------------------------------------------------
# solver.convg_pwl_fd([20,40,80,160])

