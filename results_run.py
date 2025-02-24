# -*- coding: utf-8 -*-
"""
Created on Tue Nov  5 08:26:55 2024

@author: sarah
"""
import numpy as np
import graphics

import stokes_examples as s_exs
import stokes_control as s_control


import reyn_control as r_control
import reyn_examples as r_exs

# # BFS at Re=0, Q=2, U=0, L=4, varying H/h --> resistance error, reattachment lengths
h=1
L=4
l=1
U=0
Q_stokes=2
H = np.asarray([2.75,2.5,2.25,2,1.5,1.25,1.125])
N_exs=7

BFS_stokes = np.asarray([s_exs.BFS_H2p75L4_Re0_Q2_U0,s_exs.BFS_H2p5L4_Re0_Q2_U0,s_exs.BFS_H2p25L4_Re0_Q2_U0,
            s_exs.BFS_H2L4_Re0_Q2_U0, s_exs.BFS_H1p5L4_Re0_Q2_U0, s_exs.BFS_H1p25L4_Re0_Q2_U0,
            s_exs.BFS_H1p125L4_Re0_Q2_U0])
Ns=np.asarray([160,160,160,160,160,160,200])

x_r=np.zeros(N_exs)
y_r=np.zeros(N_exs)

dP_stokes = np.zeros(N_exs)
for k in range(N_exs):
    example = BFS_stokes[k]
    N = Ns[k]
    s_solver = s_control.Stokes_Solver(example)
    dP_stokes[k] = s_solver.get_dP(N)
    x_r[k], y_r[k] = s_solver.get_attachments(N)

graphics.plot_2D_multi([y_r,x_r], H, '', ['$y_r$','$x_r$'], ['$\mathcal{H}$','length: $\Delta x = \Delta y$'], loc='lower', colors='pri')


Q_reyn = np.zeros(N_exs)
for k in range(N_exs):
    dp = dP_stokes[k]
    args = [H[k], l]
    r_solver = r_control.Reynolds_Solver(r_exs.BFS, U, dp, args)
    reyn_pres, reyn_vel = r_solver.pwl_solve(Ns[0], plot=False)
    Q_reyn[k]=reyn_vel.flux
    
res_reyn = dP_stokes/-Q_reyn
res_stokes = dP_stokes/-Q_stokes
res_err_pct = 100*np.abs(res_reyn-res_stokes)/np.abs(res_stokes)
graphics.plot_2D_multi([res_stokes,res_reyn], H, '', ['Stokes','Reyn.'],['$\mathcal{H}$','$\mathcal{R} = -\Delta \mathcal{P} / \mathcal{Q}$'], loc='upper', colors='pri')
graphics.plot_2D(res_err_pct, H, '', ['$\mathcal{H}$','$|\mathcal{R}_{Stokes}-\mathcal{R}_{Reyn.}|/\mathcal{R}_{Stokes}\%$'], color='darkmagenta')




# #------------------------------------------------------------------------------------------------------------

# dBFS at Re=0, Q=2, U=0, L=4, H=2, varying delta --> resistance error
U=0
L=4

Q_stokes = 2

H=2
deltas_H2=np.asarray([1, 0.75, 0.5, 0.25, 0.125, 0])
N_exs=6
dBFS_H2_stokes = np.asarray([s_exs.dBFS_H2L4_d1_Re0_Q2_U0,s_exs.dBFS_H2L4_d0p75_Re0_Q2_U0,s_exs.dBFS_H2L4_d0p5_Re0_Q2_U0,
                          s_exs.dBFS_H2L4_d0p25_Re0_Q2_U0, s_exs.dBFS_H2L4_d0p125_Re0_Q2_U0, s_exs.dBFS_H2L4_d0_Re0_Q2_U0])
Ns=np.asarray([160,160,160,160,200,160])

dP_H2_stokes = np.zeros(N_exs)
for k in range(N_exs):
    example = dBFS_H2_stokes[k]
    N = Ns[k]
    s_solver = s_control.Stokes_Solver(example)
    dP_H2_stokes[k] = s_solver.get_dP(N)
    
Q_H2_reyn = np.zeros(N_exs)
for k in range(N_exs):
    dp = dP_H2_stokes[k]
    delta=deltas_H2[k]
    if delta == 0:
        args = [H,L/2]
        r_solver = r_control.Reynolds_Solver(r_exs.BFS, U, dp, args)
    else:
        args=[H,delta]
        r_solver = r_control.Reynolds_Solver(r_exs.BFS_deltaSmooth, U, dp, args)
    reyn_pres, reyn_vel = r_solver.pwl_solve(200, plot=False)
    Q_H2_reyn[k] = reyn_vel.flux
res_H2_stokes = dP_H2_stokes/-Q_stokes
res_H2_reyn = dP_H2_stokes/-Q_H2_reyn

H=1.5
deltas_H1p5=np.asarray([0.75, 0.5, 0.25, 0.125, 0])
N_exs=5
dBFS_H1p5_stokes = np.asarray([s_exs.dBFS_H1p5L4_d0p75_Re0_Q2_U0,s_exs.dBFS_H1p5L4_d0p5_Re0_Q2_U0,s_exs.dBFS_H1p5L4_d0p25_Re0_Q2_U0,
                          s_exs.dBFS_H1p5L4_d0p125_Re0_Q2_U0, s_exs.dBFS_H1p5L4_d0_Re0_Q2_U0])
Ns=np.asarray([80,80,160,200,100])

dP_H1p5_stokes = np.zeros(N_exs)
for k in range(N_exs):
    example = dBFS_H1p5_stokes[k]
    N = Ns[k]
    s_solver = s_control.Stokes_Solver(example)
    dP_H1p5_stokes[k] = s_solver.get_dP(N)
    
Q_H1p5_reyn = np.zeros(N_exs)
for k in range(N_exs):
    dp = dP_H1p5_stokes[k]
    delta=deltas_H1p5[k]
    if delta == 0:
        args = [H,L/2]
        r_solver = r_control.Reynolds_Solver(r_exs.BFS, U, dp, args)
    else:
        args=[H,delta]
        r_solver = r_control.Reynolds_Solver(r_exs.BFS_deltaSmooth, U, dp, args)
    reyn_pres, reyn_vel = r_solver.pwl_solve(200, plot=False)
    Q_H1p5_reyn[k]=reyn_vel.flux
res_H1p5_stokes = dP_H1p5_stokes/-Q_stokes
res_H1p5_reyn = dP_H1p5_stokes/-Q_H1p5_reyn

H=1.25
deltas_H1p25=np.asarray([0.25, 0.125, 0.05, 0])
N_exs=4
dBFS_H1p25_stokes = np.asarray([s_exs.dBFS_H1p25L4_d0p25_Re0_Q2_U0, s_exs.dBFS_H1p25L4_d0p125_Re0_Q2_U0,
                                s_exs.dBFS_H1p25L4_d0p05_Re0_Q2_U0, s_exs.dBFS_H1p25L4_d0_Re0_Q2_U0])
Ns=np.asarray([80,200,80,80])

dP_H1p25_stokes = np.zeros(N_exs)
for k in range(N_exs):
    example = dBFS_H1p25_stokes[k]
    N = Ns[k]
    s_solver = s_control.Stokes_Solver(example)
    dP_H1p25_stokes[k] = s_solver.get_dP(N)
    
Q_H1p25_reyn = np.zeros(N_exs)
for k in range(N_exs):
    dp = dP_H1p25_stokes[k]
    args = [H, deltas_H1p25[k]]
    delta=deltas_H1p25[k]
    if delta == 0:
        args = [H,L/2]
        r_solver = r_control.Reynolds_Solver(r_exs.BFS, U, dp, args)
    else:
        args=[H,delta]
        r_solver = r_control.Reynolds_Solver(r_exs.BFS_deltaSmooth, U, dp, args)
    reyn_pres, reyn_vel= r_solver.pwl_solve(200, plot=False)
    Q_H1p25_reyn[k] = reyn_vel.flux
res_H1p25_stokes = dP_H1p25_stokes/-Q_stokes
res_H1p25_reyn = dP_H1p25_stokes/-Q_H1p25_reyn


graphics.plot_2D_multi_multi([res_H1p25_stokes,res_H1p25_reyn,res_H1p5_stokes,res_H1p5_reyn,res_H2_stokes,res_H2_reyn],
                              [deltas_H1p25,deltas_H1p25,deltas_H1p5,deltas_H1p5,deltas_H2,deltas_H2],
                              '',  ['$\mathcal{H}=1.25$ Stokes','$\mathcal{H}=1.25$ Reyn.','$\mathcal{H}=1.5$ Stokes','$\mathcal{H}=1.5$ Reyn.','$\mathcal{H}=2$ Stokes','$\mathcal{H}=2$ Reyn.'],
                              ['$\delta$','$\mathcal{R}=-\Delta \mathcal{P}/\mathcal{Q}$'], loc='upper', colors='pri')

res_H2_err_pct = 100*np.abs(res_H2_reyn-res_H2_stokes)/np.abs(res_H2_stokes)
res_H1p5_err_pct = 100*np.abs(res_H1p5_reyn-res_H1p5_stokes)/np.abs(res_H1p5_stokes)
res_H1p25_err_pct = 100*np.abs(res_H1p25_reyn-res_H1p25_stokes)/np.abs(res_H1p25_stokes)

graphics.plot_2D_multi_multi([res_H2_err_pct,res_H1p5_err_pct,res_H1p25_err_pct], [deltas_H2,deltas_H1p5,deltas_H1p25],
                              '', ['$\mathcal{H}=2$', '$\mathcal{H}=1.5$', '$\mathcal{H}=1.25$'],
                          ['$\delta$','$|\mathcal{R}_{Reyn}-\mathcal{R}_{Stokes}|/\mathcal{R}_{Stokes}\%$'], loc='upper', colors='sec')

#------------------------------------------------------------------------------------------------------------

H = np.asarray([2.75,2.5,2.25,2,1.5,1.25])
dP_bfs = dP_stokes[:-1]


cBFS_stokes = np.asarray([s_exs.BFS_H2p75L4_noEddy_Re0_Q2_U0,s_exs.BFS_H2p5L4_noEddy_Re0_Q2_U0,s_exs.BFS_H2p25L4_noEddy_Re0_Q2_U0,
                          s_exs.BFS_H2L4_noEddy_Re0_Q2_U0,s_exs.BFS_H1p5L4_noEddy_Re0_Q2_U0,s_exs.BFS_H1p25L4_noEddy_Re0_Q2_U0])
Ns=np.asarray([100,100,100,100,100,100])
N_exs=6
dP_cbfs = np.zeros(N_exs)
for k in range(N_exs):
    example = cBFS_stokes[k]
    N = Ns[k]
    s_solver = s_control.Stokes_Solver(example)
    dP_cbfs[k] = s_solver.get_dP(N)

Q=2
dP_pct_err = 100*(abs(dP_cbfs-dP_bfs))/abs(dP_bfs)
R_pct_err=dP_pct_err/Q
graphics.plot_2D(R_pct_err, H, '', ['$\mathcal{H}$','$|\mathcal{R}_{wedge}-\mathcal{R}_{BFS}|/\mathcal{R}_{BFS}\%$'], color='darkmagenta')
# 










