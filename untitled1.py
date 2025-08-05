# -*- coding: utf-8 -*-
"""
Created on Sat Aug  2 14:35:56 2025

@author: sarah
"""

# -*- coding: utf-8 -*-
"""
Created on Thu Oct 24 13:23:27 2024

@author: sarah
"""

import reyn_control
import reyn_examples

import stokes_control
import stokes_examples

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
# Reynolds 
#------------------------------------------------------------------------------
# Reyn_Example = reyn_examples.BFS
# H=1.25  #outlet height
# l=2  #inlet length
# L=4  #total length
# args = [H, l, L]
 
Reyn_Example = reyn_examples.TriSlider
h_in=1  # inlet height
h=1/4   # apex height 
h_out = 1  #oulet height
l_in = 1  # inlet length
l_out = 1  #outlet length
l_a = 1.25  # base length A  
l_b = 0.75  # base length B
args =  [h_in, h, h_out, l_in, l_a, l_b, l_out]


# Example = examples.TriCavity
# H=4
# L=1
# h=1
# l=2     # base length
# args = [H,L,h,l]

# Example = examples.Cylinder
# r=1       # radius
# drdx = 0  # radius adjustment (shallow)
# h0 = 1/2  # clearence (max height h0+r)
# l=0.5     # inlet/outlet length
# args= [ r, h0,l, drdx]

# Reyn_Example = reyn_examples.LogisticStep
# lam = -2  # slope: -lam*(H-h)/4
# H = 2   # outlet height
# h = 1   # inlet height
# L = 2   # half length
# args = [lam,H,h, L]

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

pert = reyn_solver.fd_pert_solve(N, order=4, write=write_on, plot=plots_on, scaled=scaled_on, zoom=zoom_on, uv=uv_on, inc=inc_on, get_all=True)
e2_ps, e2_us, e2_vs =  np.nan_to_num(pert.pert2_pressure.ps_2D), pert.pert2_velocity.vx, pert.pert2_velocity.vy
e4_ps, e4_us, e4_vs =  np.nan_to_num(pert.pert4_pressure.ps_2D), pert.pert4_velocity.vx, pert.pert4_velocity.vy

#------------------------------------------------------------------------------
# Stokes solver

# args = [H, h, L*2, lam]
# Stokes_Example = stokes_examples.Logistic



args =  [h_in, h_out, h, l_in, l_a, l_b, l_out]
Stokes_Example = stokes_examples.TriSlider


#------------------------------------------------------------------------------

stokes_solver = stokes_control.Stokes_Solver(Stokes_Example, args, U, Q, Re)
stokes_ps, stokes_us, stokes_vs = stokes_solver.load(N)
stokes_ps =  np.nan_to_num(stokes_ps)
if plots_on:
    stokes_solver.load_plot(N)   

#------------------------------------------------------------------------------
# Errors 
def linf(ax, ay, bx, by):
    return np.max((np.max(np.abs(ax-bx)), np.max(np.abs(ay-by))))

def l1(ax,ay,bx,by):
    return np.sum(np.abs(ax-bx)) + np.sum(np.abs(ay-by))

def l2(ax,ay,bx,by):
    return np.sum((ax-bx)**2 + (ay-by)**2) **(1/2)

l1_stokes_V = l1(stokes_us, stokes_vs, 0, 0)
linf_stokes_V = linf(stokes_us, stokes_vs, 0, 0)
l2_stokes_V = l2(stokes_us, stokes_vs, 0, 0)


l1_stokes_P = l1(stokes_ps, 0, 0, 0)
linf_stokes_P = linf(stokes_ps, 0, 0, 0)
l2_stokes_P = l2(stokes_ps, 0, 0, 0)


#------------------------------------------------------------------------------
print('Velocity l1 % errors...')
print(l1(stokes_us, stokes_vs, reyn_us, reyn_vs)/l1_stokes_V *100)
print(l1(stokes_us, stokes_vs, adj_us, adj_vs)/l1_stokes_V *100)
print(l1(stokes_us, stokes_vs, e2_us, e2_vs)/l1_stokes_V *100)
print(l1(stokes_us, stokes_vs, e4_us, e4_vs)/l1_stokes_V *100)
print('Velocity l2 % errors...')
print(l2(stokes_us, stokes_vs, reyn_us, reyn_vs)/l2_stokes_V *100)
print(l2(stokes_us, stokes_vs, adj_us, adj_vs)/l2_stokes_V *100)
print(l2(stokes_us, stokes_vs, e2_us, e2_vs)/l2_stokes_V *100)
print(l2(stokes_us, stokes_vs, e4_us, e4_vs)/l2_stokes_V *100)
print('Velocity linf % errors...')
print(linf(stokes_us, stokes_vs, reyn_us, reyn_vs)/linf_stokes_V *100)
print(linf(stokes_us, stokes_vs, adj_us, adj_vs)/linf_stokes_V *100)
print(linf(stokes_us, stokes_vs, e2_us, e2_vs)/linf_stokes_V *100)
print(linf(stokes_us, stokes_vs, e4_us, e4_vs)/linf_stokes_V *100)


print('Pressure l1 % errors...')
print(l1(stokes_ps,0, reyn_ps, 0)/l1_stokes_P *100)
print(l1(stokes_ps, 0, adj_ps, 0)/l1_stokes_P *100)
print(l1(stokes_ps, 0, e2_ps, 0)/l1_stokes_P *100)
print(l1(stokes_ps, 0, e4_ps, 0)/l1_stokes_P *100)
print('Pressure l2 % errors...')
print(l2(stokes_ps, 0, reyn_ps, 0)/l2_stokes_P *100)
print(l2(stokes_ps, 0, adj_ps, 0)/l2_stokes_P *100)
print(l2(stokes_ps, 0, e2_ps, 0)/l2_stokes_P *100)
print(l2(stokes_ps, 0, e4_ps, 0)/l2_stokes_P *100)
print('Pressure linf % errors...')
print(linf(stokes_ps, 0, reyn_ps, 0)/linf_stokes_P *100)
print(linf(stokes_ps, 0, adj_ps, 0)/linf_stokes_P *100)
print(linf(stokes_ps, 0, e2_ps, 0)/linf_stokes_P *100)
print(linf(stokes_ps, 0, e4_ps, 0)/linf_stokes_P *100)

