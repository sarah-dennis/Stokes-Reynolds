#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug  4 15:41:51 2025

@author: sarahdennis
"""
import reyn_control
import reyn_examples
import reyn_boundary as rbc
import stokes_control
import stokes_examples

import graphics

import numpy as np

#----------------
plots_on = not True # plots p(x,y) contour-mesh and (u,v) streamlines
uv_on = False   # plots u(x,y) contour-mesh and v(x,y) contour-mesh
inc_on= False   # plots u_x + v_y contour mesh and Q(x) line
zoom_on = False # plot a zoomed frame (set params in control)
scaled_on=False # plot on scaled axis (set params in control)

#------------------------------------------------------------------------------
# boundary conditions
#------------------------------------------------------------------------------

# U: velocity {u(x,y0)=U, u(x,h(x))=0}  {v(x,y0)=0, v(x,h(x))=0} 
U =0


Q_stokes=1 # |--> dP_stokes = fixed dP

Re=0

N = 160

#------------------------------------------------------------------------------
# Errors 
def linf(ax, ay, bx, by):
    return np.max((np.max(np.abs(ax-bx)), np.max(np.abs(ay-by))))

def l1(ax,ay,bx,by):
    return np.sum(np.abs(ax-bx)) + np.sum(np.abs(ay-by))


def l2(ax,ay,bx,by):
    return np.sum((ax-bx)**2 + (ay-by)**2) **(1/2)

def get_dp(ps):
    # Ny = ps.shape[0]

    j_in = int(h_in*N/2)
    j_out = int(h_out*N/2)
    dp = ps[j_out,-1]-ps[j_in,0]
    # print(dp)
    return dp
#------------------------------------------------------------------------------

Reyn_Example = reyn_examples.Logistic
Stokes_Example = stokes_examples.Logistic
#delta = -4  # slope: -lam*(H-h)/4
h_out = 2   # outlet height
h_in = 1   # inlet height
l = 4   #  length

tests = [2, 3, 4, 6, 8, 16]
test_args = [[h_out , h_in, l, lam] for lam in tests]
e2_dps = [-20.61,-23.47,-25.35,-26.14,-26.06,-24.18] #-18.96
e4_dps = [-20.65,-23.61,-25.69,-27.28,-28.74,-66.04] #-18.95
label = '$\lambda$'
#------------------------------------------------------------------------------

# Reyn_Example = reyn_examples.TriSlider
# Stokes_Example = stokes_examples.TriSlider

# h_in=1  # inlet height
# # h=1/4   # apex height 
# h_out = 1  #oulet height
# l_in = 1  # inlet length
# l_out = 1  #outlet length
# l_a = 1.25  # base length A  
# l_b = 0.75  # base length B 

# tests = [1/8, 1/4, 1/2, 3/4, 5/4, 3/2, 7/4, 2]
# test_args = [[h_in, h0, h_out, l_in, l_a, l_b, l_out] for h0 in tests]
# label = '$h_{min}$'
#------------------------------------------------------------------------------

k = 0
num_tests=len(tests)

fun_labels= ['Reyn',  '$\epsilon^2$-PLT', '$\epsilon^4$-PLT','VA-ELT', 'TG-ELT']
num_models = 5 #reyn, VA-TG adj, e2 pert, e4 pert 

l1_V_errs= np.zeros((num_models-1,num_tests))
linf_V_errs = np.zeros((num_models-1,num_tests))
l2_V_errs = np.zeros((num_models-1,num_tests))
Q_errs = np.zeros((num_models-1, num_tests))

l1_P_errs = np.zeros((num_models,num_tests))
linf_P_errs = np.zeros((num_models,num_tests))
l2_P_errs = np.zeros((num_models,num_tests))

#------------------------------------------------------------------------------

for args in test_args:
    
#------------------------------------------------------------------------------
# Stokes 
#------------------------------------------------------------------------------
    
    stokes_solver = stokes_control.Stokes_Solver(Stokes_Example, args, U, Q_stokes, Re)
    stokes_ps, stokes_us, stokes_vs = stokes_solver.load(N)
    stokes_ps = np.nan_to_num(stokes_ps)
    stokes_dp = get_dp(stokes_ps)
    if plots_on:
        stokes_solver.load_plot(N)   
    
#------------------------------------------------------------------------------
# Reynolds 
#------------------------------------------------------------------------------
    # print(stokes_dp)
    BC = rbc.Fixed(U, stokes_dp)
    reyn_solver = reyn_control.Reynolds_Solver(Reyn_Example, BC, args)
    
    reyn_P, reyn_V = reyn_solver.fd_solve(N, plot=plots_on, scaled=scaled_on, zoom=zoom_on, uv=uv_on, inc=inc_on)
    reyn_ps = np.nan_to_num(reyn_P.ps_2D)
    reyn_Q, reyn_us, reyn_vs = reyn_V.Q, reyn_V.u, reyn_V.v
    
    adj_P, adj_V= reyn_solver.fd_adj_solve(N, plot=plots_on, scaled=scaled_on, zoom=zoom_on, uv=uv_on, inc=inc_on)
    adj_ps = np.nan_to_num(adj_P.ps_2D)
    adj_Q, adj_us, adj_vs = adj_V.Q, adj_V.u, adj_V.v
    
    adj_TG_P, adj_TG_V= reyn_solver.fd_adj_TG_solve(N, plot=plots_on, scaled=scaled_on, zoom=zoom_on, uv=uv_on, inc=inc_on)
    adj_TG_ps = np.nan_to_num(adj_TG_P.ps_2D)
    adj_TG_Q, adj_TG_us, adj_TG_vs = adj_TG_V.Q, adj_TG_V.u, adj_TG_V.v
    
    BC_e2 = rbc.Fixed(U, stokes_dp)
    reyn_solver = reyn_control.Reynolds_Solver(Reyn_Example, BC_e2, args)
    e2_P, e2_V = reyn_solver.fd_pert_solve(N, order=2, plot=plots_on, scaled=scaled_on, zoom=zoom_on, uv=uv_on, inc=inc_on)
    e2_ps = np.nan_to_num(e2_P.ps_2D)
    e2_us, e2_vs, e2_Q = e2_V.u, e2_V.v, e2_V.Q
    
    BC_e4 = rbc.Fixed(U, stokes_dp)
    reyn_solver = reyn_control.Reynolds_Solver(Reyn_Example, BC_e4, args)
    e4_P, e4_V = reyn_solver.fd_pert_solve(N, order=4, plot=plots_on, scaled=scaled_on, zoom=zoom_on, uv=uv_on, inc=inc_on)
    e4_ps = np.nan_to_num(e4_P.ps_2D)
    e4_us, e4_vs, e4_Q = e4_V.u, e4_V.v, e4_V.Q

    
    #------------------------------------------------------------------------------
    l1_stokes_V = l1(stokes_us, stokes_vs, 0, 0)
    linf_stokes_V = linf(stokes_us, stokes_vs, 0, 0)
    l2_stokes_V = l2(stokes_us, stokes_vs, 0, 0)

    l1_stokes_P = l1(stokes_ps, 0, 0, 0)
    linf_stokes_P = linf(stokes_ps, 0, 0, 0)
    l2_stokes_P = l2(stokes_ps, 0, 0, 0)
    #------------------------------------------------------------------------------
    # fun_labels= ['Reyn',  '$\epsilon^2$-PLT', '$\epsilon^4$-PLT','VA-ELT', 'TG-ELT']
    test_us = [[reyn_us, reyn_vs], [e2_us,e2_vs], [e4_us,e4_vs], [adj_us, adj_vs]]
    test_ps = [reyn_ps, e2_ps, e4_ps, adj_ps, adj_TG_ps]

    test_Qs = [reyn_Q, e2_Q, e4_Q, adj_Q]
    for i in range(len(test_us)):
        l1_V_errs[i,k] = l1(stokes_us, stokes_vs, test_us[i][0], test_us[i][1])/l1_stokes_V *100
        l2_V_errs[i,k] = l2(stokes_us, stokes_vs, test_us[i][0], test_us[i][1])/l2_stokes_V *100
        linf_V_errs[i,k] = linf(stokes_us, stokes_vs, test_us[i][0], test_us[i][1])/linf_stokes_V *100
            
        Q_errs[i,k] = np.abs(Q_stokes - test_Qs[i])/np.abs(Q_stokes)*100
    for i in range(len(test_ps)):
    
        l1_P_errs[i,k] = l1(stokes_ps, 0, test_ps[i], 0)/l1_stokes_P *100
        l2_P_errs[i,k] = l2(stokes_ps, 0, test_ps[i], 0)/l2_stokes_P *100
        linf_P_errs[i,k] = linf(stokes_ps, 0, test_ps[i], 0)/linf_stokes_P *100

    
    
    k+=1
    
    
    
graphics.plot_log_multi(l1_V_errs, tests, '$L_1$ %-error Velocity, $\Delta P = \Delta P_{Stokes}$', fun_labels, [label, '$L_1$ %-error'])
graphics.plot_log_multi(l2_V_errs, tests, '$L_2$ %-error Velocity, $\Delta P = \Delta P_{Stokes}$', fun_labels, [label, '$L_2$ %-error'])
graphics.plot_log_multi(linf_V_errs, tests, '$L_\infty$ %-error Velocity, $\Delta P = \Delta P_{Stokes}$',  fun_labels,  [label, '$L_\infty$ %-error'])
graphics.plot_log_multi(Q_errs, tests, '$Q$ error %, $\Delta P = \Delta P_{Stokes}$',  fun_labels,  [label, 'Q % error'])


graphics.plot_log_multi(l1_P_errs, tests, '$L_1$ %-error Pressure, $\Delta P = \Delta P_{Stokes}$', fun_labels, [label, '$L_1$ %-error'])
graphics.plot_log_multi(l2_P_errs, tests, '$L_2$ %-error Pressure, $\Delta P = \Delta P_{Stokes}$',  fun_labels, [label, '$L_2$ %-error'])
graphics.plot_log_multi(linf_P_errs, tests, '$L_\infty$ %-error Pressure, $\Delta P = \Delta P_{Stokes}$',  fun_labels,  [label, '$L_\infty$ %-error'])






