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
U =1


Q=0


BC = rbc.Mixed(U, Q)

Re=0


N = 160 # grid size |1|= N

#------------------------------------------------------------------------------
#TODO: select example

#------------------------------------------------------------------------------
Reyn_Example = reyn_examples.TriCavity
Stokes_Example = stokes_examples.TriCavity

# H
L = 2 
# args = [H, L]

Hs = [0.35, 0.5, 0.75, 1, 2, 3, 4]
num_tests=len(Hs)
test_args = [[H, L] for H in Hs]

exstr = 'Triangular Cavity'
label = '$H$'

#------------------------------------------------------------------------------

k = 0

l1_V_errs= np.zeros(num_tests)
linf_V_errs =  np.zeros(num_tests)
l2_V_errs =  np.zeros(num_tests)

for args in test_args:
#------------------------------------------------------------------------------
# Reynolds 
#------------------------------------------------------------------------------
    reyn_solver = reyn_solvers.Reynolds_Solver(Reyn_Example, BC, args)
    reyn_solution = reyn_solver.fd_solve(N)

    reyn_ps, reyn_us, reyn_vs = np.nan_to_num(reyn_solution.pressure.ps_2D),reyn_solution.velocity.u,reyn_solution.velocity.v
    
#------------------------------------------------------------------------------
# Stokes 
#------------------------------------------------------------------------------
    
    stokes_solver = stokes_control.Stokes_Solver(Stokes_Example, args, U, Q, Re)
    stokes_ps, stokes_us, stokes_vs, stokes_dp = stokes_solver.load(N)


#------------------------------------------------------------------------------
# relative norm
    l1_stokes_V = l1(stokes_us, stokes_vs, 0, 0)
    linf_stokes_V = linf(stokes_us, stokes_vs, 0, 0)
    l2_stokes_V = l2(stokes_us, stokes_vs, 0, 0)

#------------------------------------------------------------------------------


    l1_V_errs[k] = l1(stokes_us, stokes_vs,  reyn_us, reyn_vs)/l1_stokes_V *100
    l2_V_errs[k] = l2(stokes_us, stokes_vs, reyn_us, reyn_vs)/l2_stokes_V *100
    linf_V_errs[k] = linf(stokes_us, stokes_vs,  reyn_us, reyn_vs)/linf_stokes_V *100

    k+=1
    
    
graphics.plot_2D(l2_V_errs, Hs, f'Velocity rel. %-error, {exstr}', [label, 'rel. %-error'], color='forestgreen')


