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
plots_on = False + True
uv_on =  not True # plot u(x,y) & v(x,y)
inc_on=  not True # plot ux + vy =? 0
zoom_on = not True    # plot a zoomed-in window, set location in reyn_control.py
scaled_on= False  # plot in scaled variables x/X, y/Y etc.

#------------------------------------------------------------------------------
## Piecewise-linear examples 
##       (analytic or finite difference solution)
#------------------------------------------------------------------------------

# Example = examples.BFS
# H=1
# h=2
# l=2
# L=4
# args =  [h, H, l, L]

# Example = examples.BFS_2
# H=1.5
# h=1
# l=1
# L=3
# args =  [h, H, l, L]

# Example = examples.pwl_cont_wave
# H = 2
# L = 4
# args = [H,L]


# Example = examples.BFS_deltaSmooth
# H = 2
# delta = 1
# args = [H,delta]


# Example = examples.linear 
# h0 = 1
# m = 0
# L = 4
# args = [h0,m,L]

# Example = examples.TriSlider
# h_in=1
# h=2
# h_out = h_in
# l_in = 1
# l_out = 1
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
#------------------------------------------------------------------------------
Example = examples.Sinusoid
H=1
delta = 1/4
L=4 #k=2pi/L
args = [H, delta, L]

# Example = examples.LambdaBump 
# lam=-1/2
# H=1
# l=2
# h0 = 0.5
# args=[lam, H, l, h0] 


# Example = examples.Cylinder
# r = 1      # radius
# h0 = 1/2   # cleareance: Hin=Hout = h0 + r 
# l = 0.5    # length: in let=outlet 
# drdx = 0   # depth: inlet=outlet = l6+drdx & Hin=Hout = h0+r-drdx
# args= [ r, h0,l, drdx]


# Example = examples.Logistic
# delta = 8 # max slope: delta*(H-h)/4
# H = 2       # outlet height
# h = 1       # inlet height
# L = 4       # total length
# args = [ H, h, L, delta]


#------------------------------------------------------------------------------
# boundary conditions
#------------------------------------------------------------------------------

## U: velocity BC {u(x,y0)=U, u(x,h(x))=0}  {v(x,y0)=0, v(x,h(x))=0} 
U = 1

# fixed pressure BC {p(x0,y)=-dP, p(xL,y)=0} 
# dP = 0
# BC = bc.Fixed(U,dP)

# mixed pressure BC {dp/dx (x0,y) ~ Q, p(xL,y)=0}
Q = 1
BC = bc.Mixed(U, Q)

#------------------------------------------------------------------------------

solver = control.Reynolds_Solver(Example, BC, args)

#------------------------------------------------------------------------------
# solution methods (plots  and returns pressure, velocity )


N = 100
solver.fd_solve(N, plot=plots_on, scaled=scaled_on, zoom=zoom_on, uv=uv_on, inc=inc_on)
# print(ps.ps_1D[0], ps.ps_1D[50], ps.ps_1D[100], ps.ps_1D[150],ps.ps_1D[200],ps.ps_1D[250],ps.ps_1D[300])
# solver.pwc_schur_solve(N, plot=plots_on, scaled=scaled_on, zoom=zoom_on, uv=uv_on, inc=inc_on)

solver.pwl_schur_solve(N, plot=plots_on, scaled=scaled_on, zoom=zoom_on, uv=uv_on, inc=inc_on)


# if __name__ == '__main__':
#     solver.pwc_schur_parallel_solve(N, plot=plots_on, scaled=scaled_on, zoom=zoom_on, uv=uv_on, inc=inc_on)

solver.pwl_gmres_solve(N, plot=plots_on, scaled=scaled_on, zoom=zoom_on, uv=uv_on, inc=inc_on)
# print(ps.ps_1D[0], ps.ps_1D[50], ps.ps_1D[100], ps.ps_1D[150],ps.ps_1D[200],ps.ps_1D[250],ps.ps_1D[300])
# solver.fd_adj_TG_solve(N, plot=plots_on, scaled=scaled_on, zoom=zoom_on, uv=uv_on, inc=inc_on)
# solver.fd_adj_solve(N, plot=plots_on, scaled=scaled_on, zoom=zoom_on, uv=uv_on, inc=inc_on)
# solver.fd_pert_solve(N, order=4,  plot=plots_on, scaled=scaled_on, zoom=zoom_on, uv=uv_on, inc=inc_on)

# #------------------------------------------------------------------------------
# tests = 12

# dPs_err = np.zeros(tests)
# l1Ps_err = np.zeros(tests)
# l2Ps_err = np.zeros(tests)
# linfPs_err = np.zeros(tests)

# pwc_schur_times= np.zeros(tests)
# pwl_gmres_times=np.zeros(tests)
# pwl_schur_times=np.zeros(tests)
# fd_times= np.zeros(tests)

# Ns = np.zeros(tests)
# k_0=0
# for k in range(tests):

#     N = 2**(k_0+k+1)
#     Ns[k]=N
#     print(f'k={k+1:d} of {tests:d}, N={N:d}')
#     fd_p, fd_v, fd_t = solver.fd_solve(N, plot=plots_on, scaled=scaled_on, zoom=zoom_on, uv=uv_on, inc=inc_on)
#     pwc_schur_p, pwc_schur_v, pwc_schur_t = solver.pwc_schur_solve(N, plot=plots_on, scaled=scaled_on, zoom=zoom_on, uv=uv_on, inc=inc_on)
#     # pwl_gmres_p, pwl_gmres_v, pwl_gmres_t = solver.pwl_gmres_solve(N, plot=plots_on, scaled=scaled_on, zoom=zoom_on, uv=uv_on, inc=inc_on)
#     pwl_schur_p, pwl_schur_v, pwl_schur_t = solver.pwl_schur_solve(N, plot=plots_on, scaled=scaled_on, zoom=zoom_on, uv=uv_on, inc=inc_on)

    

#     fd_times[k] = fd_t
#     # pwl_gmres_times[k]=pwl_gmres_t
#     pwl_schur_times[k] = pwl_schur_t
#     pwc_schur_times[k] = pwc_schur_t
    
#     # dPs_err[k] = abs(fd_p.dP - pwc_schur_p.dP) 
#     # l1Ps_err[k] = sum(abs(fd_p.ps_1D-pwc_schur_p.ps_1D))/N 
#     # l2Ps_err[k] = (sum((fd_p.ps_1D-pwc_schur_p.ps_1D)**2)/N)**(1/2) 
#     # linfPs_err[k] = max(abs(fd_p.ps_1D-pwc_schur_p.ps_1D))
    
#     # dPs_err[k] = abs(pwl_gmres_p.dP - fd_p.dP) 
#     # l1Ps_err[k] =sum(abs(pwl_gmres_p.ps_1D-fd_p.ps_1D))/N 
#     # l2Ps_err[k] = (sum((pwl_gmres_p.ps_1D-fd_p.ps_1D)**2)/N)**(1/2)  
#     # linfPs_err[k] = max(abs(pwl_gmres_p.ps_1D-fd_p.ps_1D)) 
    
#     # dPs_err[k] = abs(pwl_schur_p.dP - pwc_schur_p.dP) 
#     # l1Ps_err[k] =sum(abs(pwl_schur_p.ps_1D-pwc_schur_p.ps_1D))/N 
#     # l2Ps_err[k] = (sum((pwl_schur_p.ps_1D-pwc_schur_p.ps_1D)**2)/N)**(1/2)  
#     # linfPs_err[k] = max(abs(pwl_schur_p.ps_1D-pwc_schur_p.ps_1D)) 


#     # dPs_err[k] = abs(pwl_schur_p.dP - fd_p.dP) 
#     # l1Ps_err[k] =sum(abs(pwl_schur_p.ps_1D-fd_p.ps_1D))/N 
#     # l2Ps_err[k] = (sum((pwl_schur_p.ps_1D-fd_p.ps_1D)**2)/N)**(1/2)  
#     # linfPs_err[k] = max(abs(pwl_schur_p.ps_1D-fd_p.ps_1D)) 

#     dPs_err[k] = abs(pwc_schur_p.dP - pwl_schur_p.dP) 
#     l1Ps_err[k] =sum(abs(pwc_schur_p.ps_1D-pwl_schur_p.ps_1D))/N 
#     l2Ps_err[k] = (sum((pwc_schur_p.ps_1D-pwl_schur_p.ps_1D)**2)/N)**(1/2)  
#     linfPs_err[k] = max(abs(pwc_schur_p.ps_1D-pwl_schur_p.ps_1D)) 


# graphics.plot_log_multi([dPs_err, l1Ps_err, l2Ps_err, linfPs_err], Ns, 'Convergence Pressure Error', ['dP err', '$l_1$ $p$ error', '$l_2$ $p$ err', '$l_\infty$ $p$ err'], ['N', 'error'], log_x=True, loc='upper', bigO_on=True)

# # graphics.plot_2D_multi([fd_times, pwc_schur_times, pwl_gmres_times, pwl_schur_times], Ns, 'Run Time', ['FD', 'pwc schur', 'pwl gmres', 'pwl schur'], ['N', 'run time'], loc='right')

# graphics.plot_2D_multi([fd_times, pwc_schur_times, pwl_schur_times], Ns, 'Run Time', ['FD', 'pwc schur', 'pwl schur'], ['N', 'run time'], loc='left')
