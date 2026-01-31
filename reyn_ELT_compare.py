#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug  4 15:41:51 2025

@author: sarahdennis
"""
import reyn_control
import reyn_solvers
import reyn_examples
import reyn_boundary as rbc
import stokes_control
import stokes_examples

import graphics

import numpy as np
#------------------------------------------------------------------------------
# Errors 
def linf(ax, ay, bx, by):
    return np.max((np.max(np.abs(ax-bx)), np.max(np.abs(ay-by))))

def linf_(ax, ay, bx, by):
    norm_x = np.max(np.abs(ax-bx))
    norm_y = np.max(np.abs(ay-by))
    return np.max((norm_x, norm_y))


def l1(ax,ay,bx,by):
    return np.sum(np.abs(ax-bx)) + np.sum(np.abs(ay-by))
def l2(ax,ay,bx,by):
    return np.sum((ax-bx)**2 + (ay-by)**2) **(1/2)

#----------------
plots_on =  not True # plots p(x,y) contour-mesh and (u,v) streamlines
uv_on = False   # plots u(x,y) contour-mesh and v(x,y) contour-mesh
inc_on= not True   # plots u_x + v_y contour mesh and Q(x) line
zoom_on = False # plot a zoomed frame (set params in control)
scaled_on=False # plot on scaled axis (set params in control)

#------------------------------------------------------------------------------
# boundary conditions
#------------------------------------------------------------------------------

# U: velocity {u(x,y0)=U, u(x,h(x))=0}  {v(x,y0)=0, v(x,h(x))=0} 
U =0


Q=1


BC = rbc.Mixed(U, Q)

Re=0


N = 20 # grid size |1|= N

#------------------------------------------------------------------------------
#TODO: select example
#------------------------------------------------------------------------------
# Reyn_Example = reyn_examples.Logistic
# Stokes_Example= stokes_examples.Logistic
# h_in = 2   # inlet height
# h_out= 1   # outlet height
# l = 16   #  length

# tests = [2, 3, 4, 6, 8, 16, 32]
# test_args = [[h_in , h_out, l, lam] for lam in tests]
# exstr = 'Logistic Step'
# label = '$\lambda$'
#------------------------------------------------------------------------------
Reyn_Example = reyn_examples.TriSlider
Stokes_Example = stokes_examples.TriSlider

h_in=1  # inlet height
h0=1/4   # apex height 
h_out = 1  #oulet height
l_in = 7  # inlet length
l_out = 7  #outlet length
l_a = 1.25  # base length A  
l_b = 0.75  # base length B 

#test h0

tests = [1/16, 1/8, 1/4, 1/2, 3/4, 5/4, 3/2, 7/4, 2]
test_args = [[h_in, h0, h_out, l_in, l_a, l_b, l_out] for h0 in tests]

exstr = 'Triangular Slider'
label = '$H_v$'

 #------------------------------------------------------------------------------
k = 0
num_tests=len(test_args)

fun_labels= ['Reyn',  '$\epsilon^2$-PLT', '$\epsilon^4$-PLT','VA-ELT', 'TG-ELT']
num_models = 5 #reyn, VA-TG adj, e2 pert, e4 pert

l1_V_errs= np.zeros((num_models,num_tests))
linf_V_errs = np.zeros((num_models,num_tests))
l2_V_errs = np.zeros((num_models,num_tests))


l1_P_errs = np.zeros((num_models,num_tests))
linf_P_errs = np.zeros((num_models,num_tests))
l2_P_errs = np.zeros((num_models,num_tests))

for args in test_args:
#------------------------------------------------------------------------------
# Reynolds 
#------------------------------------------------------------------------------
    reyn_solver = reyn_solvers.Reynolds_Solver(Reyn_Example, BC, args)
    reyn_solution = reyn_solver.fd_solve(N)
    
    reyn_ps,reyn_us, reyn_vs = np.nan_to_num(reyn_solution.pressure.ps_2D),reyn_solution.velocity.u,reyn_solution.velocity.v
    
    VA_ELT_solution = reyn_solver.fd_VA_ELT_solve(N)
    VA_ELT_ps,VA_ELT_us, VA_ELT_vs = np.nan_to_num(VA_ELT_solution.pressure.ps_2D),VA_ELT_solution.velocity.u,VA_ELT_solution.velocity.v
    
    TG_ELT_solution = reyn_solver.fd_TG_ELT_solve(N)
    TG_ELT_ps,TG_ELT_us,TG_ELT_vs = np.nan_to_num(TG_ELT_solution.pressure.ps_2D),TG_ELT_solution.velocity.u,TG_ELT_solution.velocity.v

    pert2_solution, pert4_solution = reyn_solver.fd_pert_solve(N, order=4, get_both=True)
    e2_ps, e2_us, e2_vs =  np.nan_to_num(pert2_solution.pressure.ps_2D), pert2_solution.velocity.u, pert2_solution.velocity.v

    e4_ps, e4_us, e4_vs =  np.nan_to_num(pert4_solution.pressure.ps_2D), pert4_solution.velocity.u, pert4_solution.velocity.v

#------------------------------------------------------------------------------
# Stokes 
#------------------------------------------------------------------------------
    
    stokes_solver = stokes_control.Stokes_Solver(Stokes_Example, args, U, Q, Re)
    stokes_ps, stokes_us, stokes_vs, stokes_dp = stokes_solver.load(N)
    stokes_ps = np.nan_to_num(stokes_ps)

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
    # fun_labels= ['Reyn',  '$\epsilon^2$-PLT', '$\epsilon^4$-PLT','VA-ELT', 'TG-ELT']
    test_us = [[reyn_us, reyn_vs], [e2_us,e2_vs], [e4_us,e4_vs], [VA_ELT_us, VA_ELT_vs],  [TG_ELT_us, TG_ELT_vs]]
    test_ps = [reyn_ps, e2_ps, e4_ps, VA_ELT_ps, TG_ELT_ps]
    
    for i in range(len(test_us)):
        # l1_V_errs[i,k] = l1(stokes_us, stokes_vs, test_us[i][0], test_us[i][1])/l1_stokes_V *100
        l2_V_errs[i,k] = l2(stokes_us, stokes_vs, test_us[i][0], test_us[i][1])/l2_stokes_V *100
        # linf_V_errs[i,k] = linf(stokes_us, stokes_vs, test_us[i][0], test_us[i][1])/linf_stokes_V *100
    
    for i in range(len(test_ps)):
    
        # l1_P_errs[i,k] = l1(stokes_ps, 0, test_ps[i], 0)/l1_stokes_P *100
        l2_P_errs[i,k] = l2(stokes_ps, 0, test_ps[i], 0)/l2_stokes_P *100
        # linf_P_errs[i,k] = linf(stokes_ps, 0, test_ps[i], 0)/linf_stokes_P *100
 
    k+=1
    
    
# graphics.plot_log_multi(l1_V_errs[L:-1], tests, f'$L_1$ rel. %-error Velocity, {exstr} $Q=${Q:.1f}', fun_labels, [label, '$L_1$ rel. %-error'],loc='left')
# graphics.plot_log_multi(l2_V_errs[:-1], tests, f'Velocity $L_2$ rel. %-error, {exstr}', fun_labels, [label, 'Velocity $L_2$ rel. %-error'],loc='center')#,loc='left')
graphics.plot_2D_multi(l2_V_errs[:-1], tests, f'Velocity $L_2$ rel. %-error, {exstr}', fun_labels, [label, 'Velocity $L_2$ rel. %-error'],loc='right')#,loc='right')
# graphics.plot_log_multi(linf_V_errs, tests, f'$L_\infty$ rel. %-error Velocity, {exstr} $Q=${Q:.1f}',  fun_labels,  [label, '$L_\infty$ rel. %-error'],loc='left')


# graphics.plot_log_multi(l1_P_errs, tests, f'$L_1$ rel. %-error Pressure, {exstr} $Q=${Q:.1f}', fun_labels, [label, '$L_1$ rel. %-error'],loc='left')
# graphics.plot_log_multi(l2_P_errs, tests, f'Pressure $L_2$ rel. %-error, {exstr}',  fun_labels, [label, 'Presure $L_2$ rel. %-error '],loc='right')#,loc='left')
graphics.plot_2D_multi(l2_P_errs, tests, f'Pressure $L_2$ rel. %-error, {exstr}',  fun_labels, [label, 'Pressure $L_2$ rel. %-error '],loc='left')#,loc='lower')

# graphics.plot_log_multi(linf_P_errs, tests, f'$L_\infty$ rel. %-error Pressure, {exstr} $Q=${Q:.1f}',  fun_labels,  [label, '$L_\infty$ rel. %-error'],loc='left')


