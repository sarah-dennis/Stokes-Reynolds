#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug  4 15:41:51 2025

@author: sarahdennis
"""
import reyn_control
import reyn_examples

import stokes_control
import stokes_examples

import graphics

import numpy as np

#----------------
plots_on = True
uv_on = False
inc_on=False
zoom_on = False 
write_on = False
scaled_on=False

#------------------------------------------------------------------------------
# boundary conditions
#------------------------------------------------------------------------------

# U: velocity {u(x,y0)=U, u(x,h(x))=0}  {v(x,y0)=0, v(x,h(x))=0} 
U =0
# dP: 1D pressure {p(x0,y)=, u(x,h(x))=0} 

reyn_dP =0

Q=1

Re=0

N = 80
#------------------------------------------------------------------------------
# Errors 
def linf(ax, ay, bx, by):
    return np.max((np.max(np.abs(ax-bx)), np.max(np.abs(ay-by))))

def l1(ax,ay,bx,by):
    return np.sum(np.abs(ax-bx)) + np.sum(np.abs(ay-by))

def l2(ax,ay,bx,by):
    return np.sum((ax-bx)**2 + (ay-by)**2) **(1/2)

#------------------------------------------------------------------------------

Reyn_Example = reyn_examples.Logistic
Stokes_Example= stokes_examples.Logistic
#delta = -4  # slope: -lam*(H-h)/4
H = 2   # outlet height
h = 1   # inlet height
l = 4   #  length

tests = [-1, -2, -4, -8]

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

# tests = [1/16, 1/8, 1/4, 1/2]
#------------------------------------------------------------------------------

k = 0
num_tests=len(tests)

fun_labels= ['reyn', 'VA-TG', 'TG', '$\epsilon^2$ pert', '$\epsilon^4$ pert']
num_models = 5 #reyn, VA-TG adj, e2 pert, e4 pert 

l1_V_errs= np.zeros((num_models,num_tests))
linf_V_errs = np.zeros((num_models,num_tests))
l2_V_errs = np.zeros((num_models,num_tests))


l1_P_errs = np.zeros((num_models,num_tests))
linf_P_errs = np.zeros((num_models,num_tests))
l2_P_errs = np.zeros((num_models,num_tests))

#------------------------------------------------------------------------------

# label = '$h_{min}$'
# for h in tests:
#     args =  [h_in, h, h_out, l_in, l_a, l_b, l_out]


label = '$\delta$'
for delta in tests:
    args = [ H, h, l, delta]
#------------------------------------------------------------------------------
# Reynolds 
#------------------------------------------------------------------------------
    
    reyn_solver = reyn_control.Reynolds_Solver(Reyn_Example, U, reyn_dP, args)
    
    reyn_P, reyn_V = reyn_solver.fd_solve(N, plot=plots_on, scaled=scaled_on, zoom=zoom_on, uv=uv_on, inc=inc_on)
    reyn_ps = np.nan_to_num(reyn_P.ps_2D)
    reyn_us = reyn_V.vx
    reyn_vs = reyn_V.vy
    
    adj_P, adj_V= reyn_solver.fd_adj_solve(N, write_on, plot=plots_on, scaled=scaled_on, zoom=zoom_on, uv=uv_on, inc=inc_on, reynFlux=True)
    adj_ps = np.nan_to_num(adj_P.ps_2D)
    adj_us = adj_V.vx
    adj_vs = adj_V.vy
    
    adj_TG_P, adj_TG_V= reyn_solver.fd_adj_TG_solve(N, write_on, plot=plots_on, scaled=scaled_on, zoom=zoom_on, uv=uv_on, inc=inc_on)
    adj_TG_ps = np.nan_to_num(adj_TG_P.ps_2D)
    adj_TG_us = adj_TG_V.vx
    adj_TG_vs = adj_TG_V.vy
    
    pert = reyn_solver.fd_pert_solve(N, order=4, write=write_on, plot=plots_on, scaled=scaled_on, zoom=zoom_on, uv=uv_on, inc=inc_on, get_all=True)
    e2_ps, e2_us, e2_vs =  np.nan_to_num(pert.pert2_pressure.ps_2D), pert.pert2_velocity.vx, pert.pert2_velocity.vy
    e4_ps, e4_us, e4_vs =  np.nan_to_num(pert.pert4_pressure.ps_2D), pert.pert4_velocity.vx, pert.pert4_velocity.vy
    
#------------------------------------------------------------------------------
# Stokes 
#------------------------------------------------------------------------------
    
    stokes_solver = stokes_control.Stokes_Solver(Stokes_Example, args, U, Q, Re)
    stokes_ps, stokes_us, stokes_vs = stokes_solver.load(N)
    stokes_ps =  np.nan_to_num(stokes_ps)
    if plots_on:
        stokes_solver.load_plot(N)   
    
    
    #------------------------------------------------------------------------------
    

    l1_stokes_V = l1(stokes_us, stokes_vs, 0, 0)
    linf_stokes_V = linf(stokes_us, stokes_vs, 0, 0)
    l2_stokes_V = l2(stokes_us, stokes_vs, 0, 0)


    l1_stokes_P = l1(stokes_ps, 0, 0, 0)
    linf_stokes_P = linf(stokes_ps, 0, 0, 0)
    l2_stokes_P = l2(stokes_ps, 0, 0, 0)
    #------------------------------------------------------------------------------
    l1_V_errs[0,k] = l1(stokes_us, stokes_vs, reyn_us, reyn_vs)/l1_stokes_V *100
    l1_V_errs[1,k] = l1(stokes_us, stokes_vs, adj_us, adj_vs)/l1_stokes_V *100
    l1_V_errs[2,k] = l1(stokes_us, stokes_vs, adj_TG_us, adj_TG_vs)/l1_stokes_V *100
    l1_V_errs[3,k] = l1(stokes_us, stokes_vs, e2_us, e2_vs)/l1_stokes_V *100
    l1_V_errs[4,k] = l1(stokes_us, stokes_vs, e4_us, e4_vs)/l1_stokes_V *100
                           
    l2_V_errs[0,k] = l2(stokes_us, stokes_vs, reyn_us, reyn_vs)/l2_stokes_V *100
    l2_V_errs[1,k] = l2(stokes_us, stokes_vs, adj_us, adj_vs)/l2_stokes_V *100
    l2_V_errs[2,k] = l2(stokes_us, stokes_vs, adj_TG_us, adj_TG_vs)/l2_stokes_V *100
    l2_V_errs[3,k] = l2(stokes_us, stokes_vs, e2_us, e2_vs)/l2_stokes_V *100
    l2_V_errs[4,k] = l2(stokes_us, stokes_vs, e4_us, e4_vs)/l2_stokes_V *100
    
    linf_V_errs[0,k] = linf(stokes_us, stokes_vs, reyn_us, reyn_vs)/linf_stokes_V *100
    linf_V_errs[1,k] = linf(stokes_us, stokes_vs, adj_us, adj_vs)/linf_stokes_V *100
    linf_V_errs[2,k] = linf(stokes_us, stokes_vs, adj_TG_us, adj_TG_vs)/linf_stokes_V *100
    linf_V_errs[3,k] = linf(stokes_us, stokes_vs, e2_us, e2_vs)/linf_stokes_V *100
    linf_V_errs[4,k] = linf(stokes_us, stokes_vs, e4_us, e4_vs)/linf_stokes_V *100

    l1_P_errs[0,k] = l1(stokes_ps, 0, reyn_ps, 0)/l1_stokes_P *100
    l1_P_errs[1,k] = l1(stokes_ps, 0, adj_ps, 0)/l1_stokes_P *100
    l1_P_errs[2,k] = l1(stokes_ps, 0, adj_TG_ps, 0)/l1_stokes_P *100
    l1_P_errs[3,k] = l1(stokes_ps, 0, e2_ps, 0)/l1_stokes_P *100
    l1_P_errs[4,k] = l1(stokes_ps, 0, e4_ps, 0)/l1_stokes_P *100
    
    l2_P_errs[0,k] = l2(stokes_ps, 0, reyn_ps, 0)/l2_stokes_P *100
    l2_P_errs[1,k] = l2(stokes_ps, 0, adj_ps, 0)/l2_stokes_P *100
    l2_P_errs[2,k] = l2(stokes_ps, 0, adj_TG_ps, 0)/l2_stokes_P *100
    l2_P_errs[3,k] = l2(stokes_ps, 0, e2_ps, 0)/l2_stokes_P *100
    l2_P_errs[4,k] = l2(stokes_ps, 0, e4_ps, 0)/l2_stokes_P *100
    
    linf_P_errs[0,k] = linf(stokes_ps, 0, reyn_ps, 0)/linf_stokes_P *100
    linf_P_errs[1,k] = linf(stokes_ps, 0, adj_ps, 0)/linf_stokes_P *100
    linf_P_errs[2,k] = linf(stokes_ps, 0, adj_TG_ps, 0)/linf_stokes_P *100
    linf_P_errs[3,k] = linf(stokes_ps, 0, e2_ps, 0)/linf_stokes_P *100
    linf_P_errs[4,k] = linf(stokes_ps, 0, e4_ps, 0)/linf_stokes_P *100

    k+=1
    
    
    
graphics.plot_2D_multi(l1_V_errs, tests, 'L1 error % Velocity', fun_labels, [label, 'L1 % error'])
graphics.plot_2D_multi(l2_V_errs, tests, 'L2 error % Velocity', fun_labels, [label, 'L2 % error'])
graphics.plot_2D_multi(linf_V_errs, tests, 'Linf error % Velocity',  fun_labels,  [label, 'Linf % error'])


graphics.plot_2D_multi(l1_P_errs, tests, 'L1 error % Pressure', fun_labels, [label, 'Linf % error'])
graphics.plot_2D_multi(l2_P_errs, tests, 'L2 error % Pressure',  fun_labels, [label, 'Linf % error'])
graphics.plot_2D_multi(linf_P_errs, tests, 'Linf error % Pressure',  fun_labels,  [label, 'Linf % error'])






