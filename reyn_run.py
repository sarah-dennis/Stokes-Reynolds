# -*- coding: utf-8 -*-
"""
Created on Thu Oct 24 13:23:27 2024

@author: sarah
"""


import reyn_boundary as bc
import reyn_examples as examples
import reyn_solvers as solvers

#-------------------plotting---------------------------------------------------

plots_on = True
uv_on = False          # plot u(x,y) & v(x,y) & |(u,v)|
inc_on = False         # plot ux + vy =? 0
zoom_on = True #False        # plot a zoomed-in window, set location in reyn_solution.py
scaled_on = False      # plot in scaled variables x/X, y/Y etc.

#------------------------------------------------------------------------------
## Piecewise-linear examples 
##       (analytic or finite difference solution)
#------------------------------------------------------------------------------

# Example = examples.BFS
# H=1 
# h=2
# l=8
# l_out=8
# args =  [h, H, l, l_out]

# Example = examples.multi_step
# H = 2
# L = 4
# args = [H,L]

# Example = examples.linear
# h0 = 1
# m = 2
# L = 4
# args =  [h0, m, L]

Example = examples.BFS_deltaSmooth
H = 2
h=1
delta = 1/2
L=16
args = [H,h,L,delta]


# Example = examples.BFS_noEddy
# h = 1
# H = 2
# l = 1
# L = 4
# xr = 0.5
# yr = 0.5
# args = [h, H, l, L, xr, yr]


# Example = examples.TriSlider
# h_in=1
# h=2
# h_out = h_in
# l_in = 7
# l_out = 7
# l_a = 1.25
# l_b = 0.75
# args =  [h_in, h, h_out, l_in, l_a, l_b, l_out]



# Example = examples.TriCavity
# H=2 # apex height
# l_a = 1.25
# l_b = 0.75
# args = [H, l_a, l_b]

#------------------------------------------------------------------------------
## Smooth examples  
##      (finite difference solution only)
# ------------------------------------------------------------------------------
# Example = examples.Sinusoid
# H=1
# delta = 1/2
# k = 1 #* 2pi
# L=2
# args = [H, delta, k, L]

# Example = examples.Sinusoid_2
# H=1
# delta = 1/4
# k = 2 # period k * pi on length 2l
# l = 1 # half length of texture
# L=3 # half length total length
# args = [H, delta, k, l, L]


# Example = examples.LambdaBump 
# lam=-1/2
# H=1
# l=2
# h0 = 0.5
# args=[lam, H, l, h0] 


# Example = examples.Cylinder
# r = 1      # radius
# h0 = 1/2   # cleareance: Hin=Hout = h0 + r 
# l = 0.5    # length: inlet=outlet 
# drdx = 0   # depth: inlet=outlet = l6+drdx & Hin=Hout = h0+r-drdx
# args= [ r, h0,l, drdx]


# Example = examples.Logistic
# delta = 8 # max slope: delta*(H-h)/4
# H = 2     # outlet height
# h = 1       # inlet height
# L = 8     # total length
# args = [ H, h, L, delta]


#------------------------------------------------------------------------------
# boundary conditions
#------------------------------------------------------------------------------

## U: velocity BC {u(x,y0)=U, u(x,h(x))=0}  {v(x,y0)=0, v(x,h(x))=0} 
U = 0

#fixed pressure BC {p(x0,y)=-dP, p(xL,y)=0} 
# dP = 8
# BC = bc.Fixed(U,dP)

# mixed pressure BC {dp/dx (x0,y) ~ Q, p(xL,y)=0}
Q = 1
BC = bc.Mixed(U, Q)

#------------------------------------------------------------------------------

solver = solvers.Reynolds_Solver(Example, BC, args)

#------------------------------------------------------------------------------
# solution methods (plots  and returns pressure, velocity )


N = 200
# solution = solver.fd_solve(N)
# 
# solution = solver.pwc_schur_solve(N)

# if __name__ == '__main__':
#     solution = solver.pwc_schur_parallel_solve(N)

solution = solver.pwl_schur_solve(N)

# solution = solver.pwl_gmres_solve(N)

# solution = solver.fd_TG_ELT_solve(N)
# 
# solution = solver.fd_VA_ELT_solve(N)

# solution = solver.fd_pert_solve(N, order=2)
# solution = solver.fd_pert_solve(N, order=4, get_both=False)

if plots_on:
    solution.p_plot(scaled=scaled_on, zoom=zoom_on)
    solution.v_plot(scaled=scaled_on, zoom=zoom_on, uv=uv_on, inc=inc_on)
#------------------------------------------------------------------------------
