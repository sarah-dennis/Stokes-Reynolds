# -*- coding: utf-8 -*-
"""
Created on Thu Oct 24 13:23:27 2024

@author: sarah
"""

import reyn_control as control
import reyn_examples as examples
#------------------------------------------------------------------------------
# Instructions: 
#------------------------------------------------------------------------------
 
# I. Uncomment and select an example. 
#    class structure: examples < reyn_examples < reyn_heights < domain
#           examples sets fluid parameters dP, U, 
#           reyn_examples sets args (geometric parameters) H, L, h, r, etc. 
#           stokes_heights sets grid size N=1/dx and discretizes the height function

#       /> N, U, dP = [100, 1 , 0]
#       /> Example = examples.BFS
#       /> args = [H=2,l=4] 

# II. Initialize the solver with the selected example
#        /> solver = control.Reynolds_Solver(Example, U, dP, args)

# IV. Choose a method from the solver, provide a grid size N to run or load 
#       the solver in reyn_control has several solution, graphing and convergence methods,

#        /> solver.fd_solve(N, plot=plots_on)    # Finite difference solve

#        /> solver.pwl_solve(N, plot=plots_on)   # Analytic solve for piecewise-linear heights

#        /> solver.fd_adj_solve(N, plot=plots_on) #Finite difference solve for *adjusted Reynolds equation*
#------------------------------------------------------------------------------

plots_on=True
zoom_on=False

# Grid size
N = 100        # N = 1/dx = 1/dy
 
# surface velocity
U=0             # u(x,0) = U; u(x,h) = 0
                # v(x,0) = 0; v(x,h) = 0
                
dP = 20          # p(0,y) = -dP; p(xf,y) = 0

# space parmaters --> args
l=1         
h = .125         
H=1.5            

xr = 0.35       
yr = 0.41
delta = 0.25
r=0.125
k=3.14

#------------------------------------------------------------------------------
# Piecewise-linear examples 
#       analytic or finite difference solution
#------------------------------------------------------------------------------
Example = examples.BFS
args = [H,l]

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

#------------------------------------------------------------------------------
# Continuous (not piecewise linear) examples 
#      finite difference solution only
#------------------------------------------------------------------------------
# Example = examples.Bump # 
# args=[r,k]

# Example = examples.CircCavity
# args=[r,h,l]

# Example = examples.TriCavity
# args = [H,l]



#------------------------------------------------------------------------------

solver = control.Reynolds_Solver(Example, U, dP, args)
# solver.fd_solve(N, plot=plots_on, zoom=zoom_on)
# solver.pwl_solve(N, plot=plots_on, zoom=zoom_on)
solver.fd_adj_solve(N, plot=plots_on, zoom=zoom_on)

# solver.compare_pwl_fd([20,40,80,160,320])
# solver.load_plot(N, zoom=zoom_on)