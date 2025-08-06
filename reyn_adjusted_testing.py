# -*- coding: utf-8 -*-
"""
Created on Tue Mar 11 09:40:38 2025

@author: sarah
"""

import reyn_control as control
import reyn_examples as examples

import numpy as np
import graphics





# Example = examples.LogisticStep

# lam = -3  # slope: -lam*(H-h)/4
# H = 2   # outlet height

# h = 1   # inlet height
# l = 2   # half length

# args = [lam,H,h, l]
#------------------------------------------------------------------------------
# Convergence to analytic solution
#------------------------------------------------------------------------------
Example = examples.Cylinder

r=2
h0 = 0.05
l=0
drdx = 0

d = h0+r
print(d/r)
args= [ r, h0,l, drdx]




U=-1/2
dP=0

    
# analytic (Wannier) solution for Cylinder

s = np.sqrt(d**2 - r**2)
log_gamma = np.log((d+s)/(d-s))
B = 2*(d+s)*U/log_gamma
C = 2*(d-s)*U/log_gamma   
F = U/log_gamma    

def delta(x,y):
    return (x**2 + y**2 + s**2)/(2*y)

def p_stokes_Wa(x,y):
    d = delta(x,y)
    
    b = -B*x*(s+y)/(y**2 * (d+s)**2)
    c = -C*x*(s-y)/(y**2 * (d-s)**2)
    f = -4*F*x*s/(y* (d**2 - s**2)) 

    return b + c + f


# def p_reyn_Wa(x):
    
#     return -8*(r**2)*U*x/(2*r*h0 + x**2)**2


def pressure_cylinder(ex, N):

    ps_stokes = np.zeros((ex.Ny, ex.Nx))
    for i in range(ex.Nx):
        x = ex.xs[i]
        h = ex.hs[i]
       
        for j in range(ex.Ny):
            y = ex.ys[j]

            if y <= h and y > 0:
                ps_stokes[j,i]= p_stokes_Wa(x,y)
            else:
                ps_stokes[j,i]=None
    return ps_stokes



N = 100
ex = Example(U, dP, N, args)



solver = control.Reynolds_Solver(Example, U, dP, args)
adj_pressure, _ = solver.fd_adj_solve(N, write=False, plot=False, reynFlux=False)
adj_pressure_old = adj_pressure.ps_2D - adj_pressure.sigmas
adj_ps = adj_pressure.ps_2D/abs(U/h0)
adj_ps_old = adj_pressure_old/abs(U/h0)

reyn_pressure, _ = solver.fd_solve(N, write=False, plot=False)
reyn_ps = reyn_pressure.ps_2D/abs(U/h0)

stokes_ps = pressure_cylinder(ex, N)/abs(U/h0)

# pert_pressure = solver.fd_pert_solve(N, 4, write=False, plot=False, get_all=True)
# pert2_ps = pert_pressure.pert2_pressure.ps_2D/abs(U/h0)
# pert4_ps = pert_pressure.pert4_pressure.ps_2D/abs(U/h0)

x_start = -r
i_start = int(x_start*N) + ex.Nx//2   

x_stop = r
i_stop = int(x_stop*N) + ex.Nx//2

y_max = max(ex.hs[i_start],ex.hs[i_stop])



graphics.plot_contour_mesh(stokes_ps, ex.xs/r, ex.ys/r, 'Wannier Analytic Stokes pressure', ['p','x','y'], vmin=-3, vmax=3)
graphics.plot_contour_mesh(adj_ps_old, ex.xs/r, ex.ys/r, 'Takeuchi-Gu Adjusted-Reynolds pressure', ['p','x','y'], vmin=-3, vmax=3)
graphics.plot_contour_mesh(adj_ps, ex.xs/r, ex.ys/r, 'New Adjusted-Reynolds pressure', ['p','x','y'], vmin=-3, vmax=3)
graphics.plot_contour_mesh(reyn_ps, ex.xs/r, ex.ys/r, 'Reynolds pressure', ['p','x','y'], vmin=-3, vmax=3)
# graphics.plot_contour_mesh(pert2_ps, ex.xs/r, ex.ys/r, '2nd Perturbed Reynolds pressure', ['p','x','y'], vmin=-3, vmax=3)
# graphics.plot_contour_mesh(pert4_ps, ex.xs/r, ex.ys/r, '4th perturbed Reynolds pressure', ['p','x','y'], vmin=-3, vmax=3)



stokes_ps_nml = stokes_ps - reyn_ps
adj_ps_nml = adj_ps - reyn_ps
adj_ps_TG_nml = adj_ps_old - reyn_ps
# pert2_ps_nml = pert2_ps - reyn_ps
# pert4_ps_nml = pert4_ps - reyn_ps

test_xs = [0.2, 0.5, 0.8]

for x in test_xs:
    x_test = x*r
    i_test = int(x_test*N) + ex.Nx//2
    p_adj_test = (adj_ps_nml[1:,i_test])
    p_adj_TG_test = (adj_ps_TG_nml[1:,i_test])
    # pert2_ps_test = (pert2_ps_nml[1:,i_test])
    # pert4_ps_test = (pert4_ps_nml[1:,i_test])
    p_stokes_test = (stokes_ps_nml[1:,i_test])
    
    graphics.plot_2D_multi([p_adj_test, p_adj_TG_test,  p_stokes_test], ex.ys[1:], f'$p(x_0,y)$, $x_0/r={ex.xs[i_test]/r:.1f}$', ['V.A.','T & G', 'Stokes'], ['y/r', '$p(x_0,y) h_0/U$'],loc='lower')
    # graphics.plot_2D_multi([p_adj_test, p_adj_TG_test, pert2_ps_test, p_stokes_test], ex.ys[1:], f'$p(x_0,y)$, $x_0={ex.xs[i_test]:.1f}$', ['V.A.','T & G','p2', 'Stokes'], ['y/r', '$p(x_0,y)\frac{h_0}{U}$'],loc='lower')
    # graphics.plot_2D_multi([p_adj_test, p_adj_TG_test, pert2_ps_test, pert4_ps_test], ex.ys[1:], f'$p(x_0,y)$, $x_0={ex.xs[i_test]:.1f}$', ['V.A.','T & G','p2', 'p4'], ['y/r', '$p(x_0,y)\frac{h_0}{U}$'],loc='lower')
    # graphics.plot_2D_multi([p_adj_test, p_adj_TG_test], ex.ys[1:], f'$p(x_0,y)$, $x_0={ex.xs[i_test]:.1f}$', ['V.A.','T & G'], ['y/r', '$p(x_0,y)\frac{h_0}{U}$'],loc='lower')

    # graphics.plot_2D_multi([p_adj_test, p_adj_TG_test, p_stokes_test], ex.ys[1:], f'$p(x_0,y)$, $x_0={ex.xs[i_test]:.1f}$', ['V.A.','T & G', 'Stokes'], ['y/r', '$p(x_0,y)\frac{h_0}{U}$'],loc='lower')

