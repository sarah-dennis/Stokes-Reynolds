# -*- coding: utf-8 -*-
"""
Created on Thu Oct 24 13:23:27 2024

@author: sarah
"""

import reyn_control as control
import reyn_examples as examples

#----------------
plots_on = True
uv_on = False
inc_on=False
zoom_on = False 
write_on = False
scaled_on=False#True
#------------------------------------------------------------------------------
## Piecewise-linear examples 
##       (analytic or finite difference solution)
#------------------------------------------------------------------------------
Example = examples.BFS
H=1.25
l=2
L=4
args = [H, l, L]
 

# Example = examples.BFS_deltaSmooth
# H = 2
# delta = 1.5
# args = [H,delta]

# Example = examples.TriSlider
# h_in=1
# h=0.5
# h_out = 1
# l_in = 1
# l_out = 1
# l_a = 1
# l_b = 1
# args =  [h_in, h, h_out, l_in, l_a, l_b, l_out]


# Example = examples.TriCavity
# H=4
# L=1
# h=1
# l=2
# args = [H,L,h,l]

#------------------------------------------------------------------------------
## Smooth examples  
##      (finite difference solution only)
#------------------------------------------------------------------------------
# Example = examples.Sinusoid

# H=-0.95
# h = 1
# L=4
# args = [H,h, L]

# Example = examples.LambdaBump 

# lam=-1
# H=1/4 
# l=1
# h0 = 1
# args=[lam, H, l, h0] 



# Example = examples.Cylinder
# drdx = 0
# r=1
# h0 = 1/2
# l=0.5
# args= [ r, h0,l, drdx]


Example = examples.LogisticStep

lam = -16  # slope: -lam*(H-h)/4
H = 2   # outlet height

h = 1   # inlet height
l = 2   # half length

args = [lam,H,h, l]

#------------------------------------------------------------------------------
# boundary conditions
#------------------------------------------------------------------------------

# U: velocity {u(x,y0)=U, u(x,h(x))=0}  {v(x,y0)=0, v(x,h(x))=0} 
U =0
# dP: 1D pressure {p(x0,y)=, u(x,h(x))=0} 

dP =-26.2

#------------------------------------------------------------------------------
# solution methods (plots  and returns pressure, velocity )


N = 200


solver = control.Reynolds_Solver(Example, U, dP, args)
solver.fd_solve(N, plot=plots_on, scaled=scaled_on, zoom=zoom_on, uv=uv_on, inc=inc_on)
# solver.pwl_solve(N, plot=plots_on, scaled=scaled_on, zoom=zoom_on)

solver.fd_adj_solve(N, write_on, plot=plots_on, scaled=scaled_on, zoom=zoom_on, uv=uv_on, inc=inc_on)
solver.fd_pert_solve(N, order=4, write=write_on, plot=plots_on, scaled=scaled_on, zoom=zoom_on, uv=uv_on, inc=inc_on)

#------------------------------------------------------------------------------
# solver.load_plot(N, zoom=zoom_on)

#------------------------------------------------------------------------------
# solver.convg_pwl_fd([20,40,80,160])

