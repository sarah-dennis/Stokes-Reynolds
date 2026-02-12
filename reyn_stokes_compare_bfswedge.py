# -*- coding: utf-8 -*-
"""
Created on Fri Jan 30 15:38:10 2026

@author: sarah
"""


import numpy as np
import stokes_control
import stokes_examples
import reyn_solvers
import reyn_examples
import graphics  
import reyn_boundary as rbc


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

#------------------------------------------------------------------------------
# boundary conditions
#------------------------------------------------------------------------------

# U: velocity {u(x,y0)=U, u(x,h(x))=0}  {v(x,y0)=0, v(x,h(x))=0} 
U =0


Q=1


BC = rbc.Mixed(U, Q)

Re=0


N = 80 # grid size |1|= N

#------------------------------------------------------------------------------
#TODO: select example

#------------------------------------------------------------------------------

# xr = 0.35
# yr = 0.4
# xr_0p75 = 0.2625
# yr_0p75 = 0.3
# xr_0p5 = 0.175
# yr_0p5 = 0.2

# args = [h_in, h_out, l_in, l_out, xr, yr]
h_in = 2
h_out = 1
l_in = 8
l_out=8

xrs = [0.356, 0.2625, 0.175, 0]
xryrs = [(0.35, 0.4), (0.2625, 0.3), (0.175, 0.2), (0,0)]
num_tests=len(xryrs)
test_args = [[h_in, h_out, l_in, l_out, xr, yr] for xr, yr in xryrs]

exstr = 'Wedged BFS'
label = '$x_w$'

#------------------------------------------------------------------------------

k = 0

l1_V_errs= np.zeros(num_tests)
linf_V_errs =  np.zeros(num_tests)
l2_V_errs =  np.zeros(num_tests)

l1_P_errs = np.zeros(num_tests)
linf_P_errs = np.zeros(num_tests)
l2_P_errs =  np.zeros(num_tests)

dP_errs =  np.zeros(num_tests)

stokes_dPs = np.zeros(num_tests)

stokes_l2Ps = np.zeros(num_tests)

reyn_dPs = np.zeros(num_tests)

for args in test_args:
    if args[-2] == 0:
        Reyn_Example = reyn_examples.BFS
        Stokes_Example = stokes_examples.BFS
        args = [h_in, h_out, l_in, l_out]
    else:
        Reyn_Example = reyn_examples.BFS_wedge
        Stokes_Example = stokes_examples.BFS_wedge

#------------------------------------------------------------------------------
# Reynolds 
#------------------------------------------------------------------------------
    reyn_solver = reyn_solvers.Reynolds_Solver(Reyn_Example, BC, args)
    reyn_solution = reyn_solver.fd_solve(N)
    reyn_dp=reyn_solution.dP
    reyn_ps, reyn_us, reyn_vs = np.nan_to_num(reyn_solution.pressure.ps_2D),reyn_solution.velocity.u,reyn_solution.velocity.v
    
    reyn_dPs[k] = reyn_dp
#------------------------------------------------------------------------------
# Stokes 
#------------------------------------------------------------------------------
    
    stokes_solver = stokes_control.Stokes_Solver(Stokes_Example, args, U, Q, Re)
    stokes_ps, stokes_us, stokes_vs, stokes_dp = stokes_solver.load(N)
    stokes_ps = np.nan_to_num(stokes_ps)
    
    stokes_dPs[k] = stokes_dp
#------------------------------------------------------------------------------
# relative norm
    l1_stokes_V = l1(stokes_us, stokes_vs, 0, 0)
    linf_stokes_V = linf(stokes_us, stokes_vs, 0, 0)
    l2_stokes_V = l2(stokes_us, stokes_vs, 0, 0)

    l1_stokes_P = l1(stokes_ps, 0, 0, 0)
    linf_stokes_P = linf(stokes_ps, 0, 0, 0)
    l2_stokes_P = l2(stokes_ps, 0, 0, 0)
    stokes_l2Ps[k] = l2_stokes_P
#------------------------------------------------------------------------------


    l1_V_errs[k] = l1(stokes_us, stokes_vs,  reyn_us, reyn_vs)/l1_stokes_V *100
    l2_V_errs[k] = l2(stokes_us, stokes_vs, reyn_us, reyn_vs)/l2_stokes_V *100
    linf_V_errs[k] = linf(stokes_us, stokes_vs,  reyn_us, reyn_vs)/linf_stokes_V *100


    l1_P_errs[k] = l1(stokes_ps, 0, reyn_ps, 0)/l1_stokes_P *100
    l2_P_errs[k] = l2(stokes_ps, 0, reyn_ps, 0)/l2_stokes_P *100
    linf_P_errs[k] = linf(stokes_ps, 0, reyn_ps, 0)/linf_stokes_P *100
 
    dP_errs[k] = abs(reyn_dp-stokes_dp)/stokes_dp *100   
 
    k+=1
    
    
# graphics.plot_2D(l2_V_errs, xrs, f'Velocity rel. %-error, {exstr}', [label, '$L_2$ rel. %-error'], color='forestgreen')

# graphics.plot_2D(l2_P_errs, h_ins, f'Pressure $L_2$ rel. %-error, {exstr}', [label, 'Pressure $L_2$ rel. %-error '])
# graphics.plot_2D(dP_errs, h_ins, f'Pressure $\Delta p$ rel. %-error, {exstr}', [label, 'Pressure $\Delta p$ rel. %-error '])
# graphics.plot_2D_multi([dP_errs, l2_P_errs], xrs, f'Pressure rel. %-error, {exstr}', ['$\Delta p$', '$|p_{Reyn}-p_{Stokes}|_2$'], [label, 'rel. %-error'],loc='left')

# graphics.plot_2D(stokes_dPs, xrs, f'Stokes Pressure Drop, {exstr}', [label, '$\Delta p$'])
# graphics.plot_2D(stokes_l2Ps, xrs, f'Stokes Pressure, {exstr}', [label, '$|p|_2$'])
graphics.plot_2D_multi([100 - 100*stokes_dPs/stokes_dPs[-1],100 - 100*stokes_l2Ps/stokes_l2Ps[-1]], xrs, f'Stokes Pressure, {exstr} vs BFS', ['$\Delta p$', '$|p|_2$'], [label, '% Change'],loc='left')


# graphics.plot_2D_multi([stokes_dPs, reyn_dPs], xrs, f'Pressure Drop, {exstr}', ['Stokes', 'Reynolds'], [label, '$\Delta p$'],loc='left')

