# -*- coding: utf-8 -*-
"""
Created on Thu Oct 24 13:23:27 2024

@author: sarah
"""

import reyn_control as control
import reyn_boundary as bc
import reyn_examples as examples

#----------------
plots_on = True
uv_on = False
inc_on= True #False
zoom_on = False 
write_on = False
scaled_on=False#True
#------------------------------------------------------------------------------
## Piecewise-linear examples 
##       (analytic or finite difference solution)
#------------------------------------------------------------------------------
# Example = examples.BFS
# H=1.25

# h=1
# l=2
# L=4
# args =  [h, H, l, L]

 

# Example = examples.BFS_deltaSmooth
# H = 2
# delta = 1.5
# args = [H,delta]


Example = examples.TriSlider
h_in=1
h=1/4
h_out = 1
l_in = 1
l_out = 1
l_a = 1.25
l_b = 0.75
args =  [h_in, h, h_out, l_in, l_a, l_b, l_out]



# Example = examples.TriCavity
# H=0.25 # apex height
# l_a = 1.25
# l_b = 0.75
# args = [H,l_a, l_b]

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

# r = 1      # radius
# h0 = 1/2   # cleareance: Hin=Hout = h0 + r 
# l = 0.5    # length: in let=outlet 
# drdx = 0   # depth: inlet=outlet = l+drdx & Hin=Hout = h0+r-drdx

# args= [ r, h0,l, drdx]


# Example = examples.Logistic

# delta = 6  # max slope: -delta*(H-h)/4
# H = 2       # outlet height
# h = 1       # inlet height
# L = 4       # total length

# args = [ H, h, L, delta]


#------------------------------------------------------------------------------
# boundary conditions
#------------------------------------------------------------------------------

## U: velocity BC {u(x,y0)=U, u(x,h(x))=0}  {v(x,y0)=0, v(x,h(x))=0} 
U = 0

## fixed pressure BC {p(x0,y)=-dP, p(xL,y)=0} 
# dP = -350#-29.61
# BC = bc.Fixed(U,dP)

## mixed pressure BC {dp/dx (x0,y) ~ Q, p(xL,y)=0}
Q = 1
BC = bc.Mixed(U, Q)

#------------------------------------------------------------------------------
# solution methods (plots  and returns pressure, velocity )


N = 320


solver = control.Reynolds_Solver(Example, BC, args)
# solver.fd_solve(N, plot=plots_on, scaled=scaled_on, zoom=zoom_on, uv=uv_on, inc=inc_on)
# solver.pwl_solve(N, plot=plots_on, scaled=scaled_on, zoom=zoom_on)
# solver.fd_adj_TG_solve(N, plot=plots_on, scaled=scaled_on, zoom=zoom_on, uv=uv_on, inc=inc_on)
solver.fd_adj_solve(N, plot=plots_on, scaled=scaled_on, zoom=zoom_on, uv=uv_on, inc=inc_on)
# solver.fd_pert_solve(N, order=4,  plot=plots_on, scaled=scaled_on, zoom=zoom_on, uv=uv_on, inc=inc_on)

#------------------------------------------------------------------------------
# solver.load_plot(N, zoom=zoom_on)

#------------------------------------------------------------------------------
# solver.convg_pwl_fd([20,40,80,160])

