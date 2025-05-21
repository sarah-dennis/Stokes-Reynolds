# -*- coding: utf-8 -*-
"""
Created on Tue Mar 11 09:40:38 2025

@author: sarah
"""


import reyn_control as control
import reyn_examples as examples

import numpy as np
import graphics


#------------------------------------------------------------------------------
# Convergence to analytic solution
#------------------------------------------------------------------------------
# Example = examples.Cylinder
# r=1/2
# h0 = 1/4
# l=1/2
# drdx = 1/8
# d = h0+r

# args= [ r, h0,l, drdx]



Example = examples.Cylinder

r=1
h0 = 0.01
l=1
d = h0+r

args= [ r, h0,l, 0]
print(f'd/a: {d/r:.2f}')
U=-1/2
dP=0

    
# analytic (Wannier) solution for Cylinder
#height = solver.example
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


def p_reyn_Wa(x):
    
    return -8*(r**2)*U*x/(2*r*h0 + x**2)**2


def pressure_cylinder(ex, N):

    ps_stokes = np.zeros((ex.Ny, ex.Nx))
    ps_reyn = np.zeros((ex.Ny, ex.Nx))    
    for i in range(ex.Nx):
        x = ex.xs[i]
        h = ex.hs[i]
        ps_reyn[:,i] = p_reyn_Wa(x) 
       
        for j in range(ex.Ny):
            y = ex.ys[j]

            if y <= h and y > 0:
                ps_stokes[j,i]= p_stokes_Wa(x,y)
            else:
                ps_stokes[j,i]=None
    return ps_stokes, ps_reyn



N = 200
ex = Example(U, dP, N, args)



solver = control.Reynolds_Solver(Example, U, dP, args)
adj_pressure, _ = solver.fd_adj_solve(N, write=False, plot=False)
reyn_pressure, _ = solver.fd_solve(N, write=False, plot=False)

adj_ps = adj_pressure.ps_2D/abs((U/h0))
reyn_ps = reyn_pressure.ps_2D/abs((U/h0))

stokes_ps,_ = pressure_cylinder(ex, N)
stokes_ps = stokes_ps/abs((U/h0))


x_start = -1
i_start = int(x_start*N) + ex.Nx//2   

x_stop = 1
i_stop = int(x_stop*N) + ex.Nx//2

y_max = max(ex.hs[i_start],ex.hs[i_stop])



graphics.plot_contour_mesh(stokes_ps, ex.xs/r, ex.ys/r, 'Wannier Analytic Stokes pressure', ['p','x','y'], vmin=-3, vmax=3)
graphics.plot_contour_mesh(adj_ps, ex.xs/r, ex.ys/r, 'Takeuchi-Gu Adjusted Reynolds pressure', ['p','x','y'], vmin=-3, vmax=3)
graphics.plot_contour_mesh(reyn_ps, ex.xs/r, ex.ys/r, 'Reynolds pressure', ['p','x','y'], vmin=-3, vmax=3)



stokes_adj = stokes_ps - reyn_ps
reyn_adj = adj_ps - reyn_ps

x_test = 0.2/d
i_test = int(x_test*N) + ex.Nx//2
p_adj_test = (reyn_adj[1:,i_test])
p_anylt_test = (stokes_adj[1:,i_test])

graphics.plot_2D_multi([p_adj_test, p_anylt_test], ex.ys[1:], f'$p(x_0,y)$, $x_0={ex.xs[i_test]:.1f}$', ['adjusted', 'stokes'], ['y', '$p(x_0,y)$'],loc='lower')

x_test = 0.5/d
i_test = int(x_test*N) + ex.Nx//2
p_adj_test = (reyn_adj[1:,i_test])
p_anylt_test = (stokes_adj[1:,i_test])

graphics.plot_2D_multi([p_adj_test, p_anylt_test], ex.ys[1:], f'$p(x_0,y)$, $x_0={ex.xs[i_test]:.1f}$', ['adjusted', 'stokes'], ['y', '$p(x_0,y)$'],loc='lower')


x_test = 0.8/d
i_test = int(x_test*N) + ex.Nx//2
p_adj_test = (reyn_adj[1:,i_test])
p_anylt_test = (stokes_adj[1:,i_test])

graphics.plot_2D_multi([p_adj_test, p_anylt_test], ex.ys[1:], f'$p(x_0,y)$, $x_0={ex.xs[i_test]:.1f}$', ['adjusted', 'stokes'], ['y', '$p(x_0,y)$'],loc='lower')




