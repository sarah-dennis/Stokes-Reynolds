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

#-------------------plotting---------------------------------------------------

plots_on = False + True
uv_on =  not True # plot u(x,y) & v(x,y)
inc_on=  not True # plot ux + vy =? 0
zoom_on =   True    # plot a zoomed-in window, set location in reyn_control.py
scaled_on= False  # plot in scaled variables x/X, y/Y etc.

#------------------------------------------------------------------------------
## Piecewise-linear examples 
##       (analytic or finite difference solution)
#------------------------------------------------------------------------------

Example = examples.BFS
H=1
h=2
l=8
L=16
args =  [h, H, l, L]

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
# L=16
# args = [H,delta,L]


# Example = examples.linear 
# h0 = 1
# m = 0
# L = 4
# args = [h0,m,L]

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
# l = 0.5    # length: in let=outlet 
# drdx = 0   # depth: inlet=outlet = l6+drdx & Hin=Hout = h0+r-drdx
# args= [ r, h0,l, drdx]


# Example = examples.Logistic
# delta = 8 # max slope: delta*(H-h)/4
# H = 2     # outlet height
# h = 1       # inlet height
# L = 16     # total length
# args = [ H, h, L, delta]


#------------------------------------------------------------------------------
# boundary conditions
#------------------------------------------------------------------------------

## U: velocity BC {u(x,y0)=U, u(x,h(x))=0}  {v(x,y0)=0, v(x,h(x))=0} 
U = 0

# fixed pressure BC {p(x0,y)=-dP, p(xL,y)=0} 
# dP = 0
# BC = bc.Fixed(U,dP)

# mixed pressure BC {dp/dx (x0,y) ~ Q, p(xL,y)=0}
Q = 1

#sinuosoid Q for DP=0
# Q = (U*H/2) * (1-(delta**2))/(1+(delta**2)/2)


BC = bc.Mixed(U, Q)

#------------------------------------------------------------------------------

solver = control.Reynolds_Solver(Example, BC, args)

#------------------------------------------------------------------------------
# solution methods (plots  and returns pressure, velocity )


N = 80
solver.fd_solve(N, plot=plots_on, scaled=scaled_on, zoom=zoom_on, uv=uv_on, inc=inc_on)

# solver.pwc_schur_solve(N, plot=plots_on, scaled=scaled_on, zoom=zoom_on, uv=uv_on, inc=inc_on)
# 
# solver.pwl_schur_solve(N, plot=plots_on, scaled=scaled_on, zoom=zoom_on, uv=uv_on, inc=inc_on)

# if __name__ == '__main__':
#     solver.pwc_schur_parallel_solve(N, plot=plots_on, scaled=scaled_on, zoom=zoom_on, uv=uv_on, inc=inc_on)

# solver.pwl_gmres_solve(N, plot=plots_on, scaled=scaled_on, zoom=zoom_on, uv=uv_on, inc=inc_on)


solver.fd_adj_TG_solve(N, plot=plots_on, scaled=scaled_on, zoom=zoom_on, uv=uv_on, inc=inc_on)

solver.fd_adj_solve(N, plot=plots_on, scaled=scaled_on, zoom=zoom_on, uv=uv_on, inc=inc_on)

solver.fd_pert_solve(N, order=4,  plot=plots_on, scaled=scaled_on, zoom=zoom_on, uv=uv_on, inc=inc_on)

#------------------------------------------------------------------------------
# tests = 12                                                                                                                                                                                                                                                              

# dPs_err = np.zeros(tests)
# l1Ps_err = np.zeros(tests)
# l2Ps_err = np.zeros(tests)
# linfPs_err = np.zeros(tests)

# pwc_schur_times= np.zeros(tests)
# pwl_schur_times=np.zeros(tests)
# fd_times= np.zeros(tests)

# dP_err_fd= np.zeros(tests)
# dP_err_pwc= np.zeros(tests)
# dP_err_pwl= np.zeros(tests)

# l2_err_fd= np.zeros(tests)
# l2_err_pwc= np.zeros(tests)
# l2_err_pwl= np.zeros(tests)

# Ns = np.zeros(tests)
# k_0=0 # start with N = 2**(k0 + 1)
# for k in range(tests):

#     N = 2**(k_0+k+1)
#     Ns[k]=N
#     print(f'k={k+1:d} of {tests:d}, N={N:d}')
#     fd_p, fd_v, fd_t = solver.fd_solve(N, plot=plots_on, scaled=scaled_on, zoom=zoom_on, uv=uv_on, inc=inc_on)
#     pwc_schur_p, pwc_schur_v, pwc_schur_t = solver.pwc_schur_solve(N, plot=plots_on, scaled=scaled_on, zoom=zoom_on, uv=uv_on, inc=inc_on)
#     pwl_schur_p, pwl_schur_v, pwl_schur_t = solver.pwl_schur_solve(N, plot=plots_on, scaled=scaled_on, zoom=zoom_on, uv=uv_on, inc=inc_on)

#     sinus_ps = solver.sinusoid_exact_sol(N, plot=plots_on, scaled=scaled_on, zoom=zoom_on, uv=uv_on, inc=inc_on)
    

    # fd_times[k] = fd_t    
    # pwc_schur_times[k] = pwc_schur_t
    # pwl_schur_times[k] = pwl_schur_t

    
    # dP_err_pwl[k] = abs(pwl_schur_p.dP) 
    # dP_err_pwc[k] = abs(pwc_schur_p.dP)
    # dP_err_fd[k] = abs(fd_p.dP)
    
    # l2_err_fd[k] = (sum((fd_p.ps_1D-sinus_ps)**2)/N)**(1/2)  
    # l2_err_pwc[k] = (sum((pwc_schur_p.ps_1D-sinus_ps)**2)/N)**(1/2)  
    # l2_err_pwl[k] = (sum((pwl_schur_p.ps_1D-sinus_ps)**2)/N)**(1/2)  

    
    # dPs_err[k] = abs(fd_p.dP - pwc_schur_p.dP) 
    # l1Ps_err[k] = sum(abs(fd_p.ps_1D-pwc_schur_p.ps_1D))/N 
    # l2Ps_err[k] = (sum((fd_p.ps_1D-pwc_schur_p.ps_1D)**2)/N)**(1/2) 
    # linfPs_err[k] = max(abs(fd_p.ps_1D-pwc_schur_p.ps_1D))

    
    # dPs_err[k] = abs(pwl_schur_p.dP) 
    # l1Ps_err[k] =sum(abs(pwl_schur_p.ps_1D-sinus_ps))/N 
    # l2Ps_err[k] = (sum((pwl_schur_p.ps_1D-sinus_ps)**2)/N)**(1/2)  
    # linfPs_err[k] = max(abs(pwl_schur_p.ps_1D-sinus_ps)) 


    # dPs_err[k] = abs(pwl_schur_p.dP - fd_p.dP) 
    # l1Ps_err[k] =sum(abs(pwl_schur_p.ps_1D-fd_p.ps_1D))/N 
    # l2Ps_err[k] = (sum((pwl_schur_p.ps_1D-fd_p.ps_1D)**2)/N)**(1/2)  
    # linfPs_err[k] = max(abs(pwl_schur_p.ps_1D-fd_p.ps_1D)) 

    # dPs_err[k] = abs(pwc_schur_p.dP - pwl_schur_p.dP) 
    # l1Ps_err[k] =sum(abs(pwc_schur_p.ps_1D-pwl_schur_p.ps_1D))/N 
    # l2Ps_err[k] = (sum((pwc_schur_p.ps_1D-pwl_schur_p.ps_1D)**2)/N)**(1/2)  
    # linfPs_err[k] = max(abs(pwc_schur_p.ps_1D-pwl_schur_p.ps_1D)) 


# graphics.plot_log_multi([dPs_err, l1Ps_err, l2Ps_err, linfPs_err], Ns, 'Convergence Pressure Error', ['dP err', '$l_1$ $p$ error', '$l_2$ $p$ err', '$l_\infty$ $p$ err'], ['N', 'error'], log_x=True, loc='upper', bigO_on=True)

# graphics.plot_2D_multi([fd_times, pwc_schur_times, pwl_schur_times], Ns, 'Run Time', ['FD', 'PWC', 'PWL'], ['$1/\Delta x$', 'run time (s)'], loc='left')

# graphics.plot_2D(pwl_schur_times, Ns, 'Run Time for PWL', ['$1/\Delta x$', 'run time (s)'], color='forestgreen', marker='s')

# graphics.plot_log_multi([dP_err_fd, dP_err_pwc, dP_err_pwl], Ns, 'Convergence in $\Delta P$: Absolute error', ['FD', 'PWC', 'PWL'], ['$1/\Delta x$', 'error'], log_x=True, loc='upper', bigO_on=True)

# graphics.plot_log_multi([l2_err_fd, l2_err_pwc, l2_err_pwl], Ns, 'Convergence in $p(x)$: $l_2$ error', ['FD', 'PWC', 'PWL'], ['$1/\Delta x$', 'error'], log_x=True, loc='upper', bigO_on=True)