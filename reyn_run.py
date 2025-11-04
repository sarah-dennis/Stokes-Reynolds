# -*- coding: utf-8 -*-
"""
Created on Thu Oct 24 13:23:27 2024

@author: sarah
"""

import reyn_control as control
import reyn_boundary as bc
import reyn_examples as examples
import graphics
import numpy as np
import time
#-------------------plotting---------------------------------------------------
plots_on =  False
uv_on =  not True # plot u(x,y) & v(x,y)
inc_on=  not True # plot ux + vy =? 0
zoom_on = not True    # plot a zoomed-in window, set location in reyn_control.py
scaled_on= False  # plot in scaled variables x/X, y/Y etc.

#------------------------------------------------------------------------------
## Piecewise-linear examples 
##       (analytic or finite difference solution)
#------------------------------------------------------------------------------

# Example = examples.BFS_2
# H=1
# h=3
# l=1
# L=3
# args =  [h, H, l, L]

 

# Example = examples.BFS_deltaSmooth
# H = 2
# delta = 1.5
# args = [H,delta]


Example = examples.TriSlider
h_in=1
h=2
h_out = h_in
l_in = 1
l_out = 1
l_a = 1.25
l_b = 0.75
args =  [h_in, h, h_out, l_in, l_a, l_b, l_out]



# Example = examples.TriCavity
# H=2 # apex height
# l_a = 1.25
# l_b = 0.75
# args = [H, l_a, l_b]

#------------------------------------------------------------------------------
## Smooth examples  
##      (finite difference solution only)
#------------------------------------------------------------------------------
# Example = examples.Sinusoid
# H=-0.95
# h = 1
# L=4
# args = [H,h, L]

Example = examples.LambdaBump 
lam=1/2
H=1
l=2
h0 = 1
args=[lam, H, l, h0] 


# Example = examples.Cylinder
# r = 1      # radius
# h0 = 1/2   # cleareance: Hin=Hout = h0 + r 
# l = 0.5    # length: in let=outlet 
# drdx = 0   # depth: inlet=outlet = l6+drdx & Hin=Hout = h0+r-drdx
# args= [ r, h0,l, drdx]


# Example = examples.Logistic
# delta = 8 # max slope: -delta*(H-h)/4
# H = 2       # outlet height
# h = 1       # inlet height
# L = 16       # total length
# args = [ H, h, L, delta]


#------------------------------------------------------------------------------
# boundary conditions
#------------------------------------------------------------------------------

## U: velocity BC {u(x,y0)=U, u(x,h(x))=0}  {v(x,y0)=0, v(x,h(x))=0} 
U = 0

# fixed pressure BC {p(x0,y)=-dP, p(xL,y)=0} 
# dP = 5
# BC = bc.Fixed(U,dP)

# # mixed pressure BC {dp/dx (x0,y) ~ Q, p(xL,y)=0}
Q = 1/2
BC = bc.Mixed(U, Q)

#------------------------------------------------------------------------------
# solution methods (plots  and returns pressure, velocity )
solver = control.Reynolds_Solver(Example, BC, args)
tests = 9
reyn_dPs = np.zeros(tests)
reyn_times = np.zeros(tests)
schur_dPs = np.zeros(tests)
schur_times = np.zeros(tests)
Ns = np.zeros(tests)
for k in range(tests):
    N = 2**(k+1)
    
    fd_p, fd_v = solver.fd_solve(N, plot=plots_on, scaled=scaled_on, zoom=zoom_on, uv=uv_on, inc=inc_on)
    reyn_dPs[k] = np.abs(fd_p.dP)
    reyn_times[k] = fd_p.time
    
    schur_p, schur_v = solver.pwc_schur_solve(N, plot=plots_on, scaled=scaled_on, zoom=zoom_on, uv=uv_on, inc=inc_on)
    schur_dPs[k] = np.abs(schur_p.dP)
    schur_times[k] = schur_p.time
    Ns[k]=N
fig_dp= graphics.plot_log_multi([reyn_dPs,schur_dPs], Ns, 'dP convergence', ['fd', 'schur'], ['N', 'dP'], loc='upper', log_x=True, log_y=True)
fig_dp.savefig('dP_fig')
fig_times = graphics.plot_log_multi([reyn_times,schur_times], Ns, 'solve time', ['fd', 'schur'], ['N', 'time'], loc='lower', log_x=True, log_y=True)
fig_times.savefig('times_fig')
# solver.pwl_gmres_solve(N, plot=plots_on, scaled=scaled_on, zoom=zoom_on, uv=uv_on, inc=inc_on)
# 
# solver.fd_adj_TG_solve(N, plot=plots_on, scaled=scaled_on, zoom=zoom_on, uv=uv_on, inc=inc_on)
# solver.fd_adj_solve(N, plot=plots_on, scaled=scaled_on, zoom=zoom_on, uv=uv_on, inc=inc_on)
# solver.fd_pert_solve(N, order=4,  plot=plots_on, scaled=scaled_on, zoom=zoom_on, uv=uv_on, inc=inc_on)
