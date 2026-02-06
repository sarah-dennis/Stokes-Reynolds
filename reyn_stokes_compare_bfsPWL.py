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

# h_in=1  # inlet height
h_out=1   # outlet height 
L = 16 # total length

h_ins = [2.75, 2, 1.25]
deltas = [1, 0.5, 0.25, 0.125, 0]
 
num_h = len(h_ins)
num_d = len(deltas)
num_tests =num_d*num_h

exstr = 'Wedge slider'
label = '$\delta$'
#------------------------------------------------------------------------------

k = 0

l1_V_errs= np.zeros((num_h,num_d))
linf_V_errs =  np.zeros((num_h,num_d))
l2_V_errs =  np.zeros((num_h,num_d))

l1_P_errs = np.zeros((num_h,num_d))
linf_P_errs = np.zeros((num_h,num_d))
l2_P_errs =  np.zeros((num_h,num_d))

dP_errs =  np.zeros((num_h,num_d))

dPs =  np.zeros((2*num_h,num_d))


for h in range(num_h):
    for d in range(num_d):
        h_in = h_ins[h]
        delta = deltas[d]
        if delta == 0:
            Reyn_Example = reyn_examples.BFS
            Stokes_Example = stokes_examples.BFS 
            
            args = [h_in, h_out, L//2, L//2]
        else:
            
            Reyn_Example = reyn_examples.BFS_pwl
            Stokes_Example = stokes_examples.BFS_pwl
            args = [h_in, h_out, L, delta]
#------------------------------------------------------------------------------
# Reynolds 
#------------------------------------------------------------------------------
        reyn_solver = reyn_solvers.Reynolds_Solver(Reyn_Example, BC, args)
        reyn_solution = reyn_solver.fd_solve(N)
        reyn_dp=reyn_solution.dP
        reyn_ps, reyn_us, reyn_vs = np.nan_to_num(reyn_solution.pressure.ps_2D),reyn_solution.velocity.u,reyn_solution.velocity.v
        dPs[h, d] = reyn_dp
#------------------------------------------------------------------------------
# Stokes 
#------------------------------------------------------------------------------
        
        stokes_solver = stokes_control.Stokes_Solver(Stokes_Example, args, U, Q, Re)
        stokes_ps, stokes_us, stokes_vs, stokes_dp = stokes_solver.load(N)
        stokes_ps = np.nan_to_num(stokes_ps)
        dPs[num_h + h, d] =stokes_dp

#------------------------------------------------------------------------------
# relative norm

        l2_stokes_V = l2(stokes_us, stokes_vs, 0, 0)
        l2_stokes_P = l2(stokes_ps, 0, 0, 0)
        
#------------------------------------------------------------------------------
        l2_V_errs[h,d] = l2(stokes_vs, 0, reyn_vs, 0)/l2_stokes_V *100
        l2_P_errs[h,d] = l2(stokes_ps, 0, reyn_ps, 0)/l2_stokes_P *100
        dP_errs[h,d] = abs(reyn_dp-stokes_dp)/stokes_dp *100   

fun_labels = ["$\mathcal{H}$ = %.2f"%h_in for h_in in h_ins] 

fun_labels_2 = [label + " Reyn" for label in fun_labels] + [label + " Stokes" for label in fun_labels] 

graphics.plot_2D_multi(l2_V_errs, deltas, f'Velocity $L_2$ rel. %-error, {exstr}',fun_labels,  [label, 'Velocity $L_2$ rel. %-error'])#,loc='right')
graphics.plot_2D_multi(l2_P_errs, deltas, f'Pressure $L_2$ rel. %-error, {exstr}', fun_labels, [label, 'Pressure $L_2$ rel. %-error '])#,loc='lower')
graphics.plot_2D_multi(dP_errs, deltas, f'Pressure $\Delta p$ rel. %-error, {exstr}',fun_labels, [label, 'Pressure $\Delta p$ rel. %-error '])#,loc='lower')
graphics.plot_2D_multi(dPs, deltas, f'Pressure $\Delta p$, {exstr}',fun_labels_2, [label, 'Pressure $\Delta p$'])#,loc='lower')

