# -*- coding: utf-8 -*-
"""
Created on Tue Nov  5 08:26:55 2024

@author: sarah
"""
import numpy as np
import graphics

# BFS at Re=0, Q=2, U=0, L=4, varying H/h --> resistance error, reattachment lengths
h=1
L=4
l=1
H = np.asarray([2.75,2.5,2.25,2,1.5,1.25,1.125])

dP_stokes = np.asarray([40.12,40.66,41.55,43.12,51.65,63.95,75.91])
dP_reyn = np.asarray([27.50,28.55,30.25,33.05,45.25,61.00,74.50])
flux = np.asarray([2,2,2,2,2,2,2])

res_stokes = dP_stokes/flux
res_reyn = dP_reyn/flux
res_err_pct = 100*np.abs(res_reyn-res_stokes)/np.abs(res_stokes)

graphics.plot_2D_multi([res_reyn, res_stokes], H, '', ['Reyn.', 'Stokes'],['$H/h$','$\mathcal{R} = -\Delta \mathcal{P} / \mathcal{Q}$'], loc='upper')
graphics.plot_2D(res_err_pct, H, '', ['$H/h$','$|\mathcal{R}_{S}-\mathcal{R}_{R}|/\mathcal{R}_{S}\%$'], color='darkmagenta')

x_r=np.asarray([0.47,0.43,0.40,0.35,0.25,0.15,0.09])
y_r=np.asarray([0.51,0.48,0.45,0.41,0.3,0.17,0.1])
graphics.plot_2D_multi([y_r,x_r], H, '', ['$y_r$','$x_r$'], ['$H/h$','$x$ , $y$'], loc='lower')

# eddy_area=0.5*x_r*y_r
# total_area= (h*l) + (H*(L-l))
# eddy_pct_area=100*(1-(total_area-eddy_area)/total_area)

# graphics.plot_2D(eddy_pct_area, H, '', ['$H/h$','%-area'])


#------------------------------------------------------------------------------------------------------------

# dBFS at Re=0, Q=2, U=0, L=4, H=2, varying delta --> resistance error

deltas_H2=np.asarray([1, 0.75, 0.5, 0.25, 0.125, 0])
deltas_H1p5=np.asarray([0.75, 0.5, 0.25, 0.125, 0])
deltas_H1p25=np.asarray([0.25, 0.125, 0.05, 0])

flux = 2


dP_H2_stokes=np.asarray([47.77, 50.85, 54.54, 59.17, 61.85, 64.19])
dP_H2_reyn=np.asarray([45.0, 47.25, 49.5, 51.75, 53.0, 54.0])

dP_H1p5_stokes=np.asarray([59.89, 61.89, 64.52, 66.64, 68.41])
dP_H1p5_reyn=np.asarray([59.0, 60.0, 61.0, 61.75, 62.25])

dP_H1p25_stokes=np.asarray([73.15, 74.24, 74.91, 75.52])
dP_H1p25_reyn=np.asarray([72.0,72.5, 72.5, 72.75])

res_H2_stokes = dP_H2_stokes/flux
res_H2_reyn = dP_H2_reyn/flux
res_H1p5_stokes = dP_H1p5_stokes/flux
res_H1p5_reyn = dP_H1p5_reyn/flux
res_H1p25_stokes = dP_H1p25_stokes/flux
res_H1p25_reyn = dP_H1p25_reyn/flux

graphics.plot_2D_multi_multi([res_H1p25_stokes,res_H1p25_reyn,res_H1p5_stokes,res_H1p5_reyn,res_H2_stokes,res_H2_reyn],[deltas_H1p25,deltas_H1p25,deltas_H1p5,deltas_H1p5,deltas_H2,deltas_H2],
                              '',  ['$H=1.25$ Stokes','$H=1.25$ Reyn.','$H=1.5$ Stokes','$H=1.5$ Reyn.','$H=2$ Stokes','$H=2$ Reyn.'],
                              ['$\delta$','$\mathcal{R}=-\Delta \mathcal{P}/\mathcal{Q}$'], loc='upper')

res_H2_err_pct = 100*np.abs(res_H2_reyn-res_H2_stokes)/np.abs(res_H2_stokes)
res_H1p5_err_pct = 100*np.abs(res_H1p5_reyn-res_H1p5_stokes)/np.abs(res_H1p5_stokes)
res_H1p25_err_pct = 100*np.abs(res_H1p25_reyn-res_H1p25_stokes)/np.abs(res_H1p25_stokes)

graphics.plot_2D_multi_multi([res_H2_err_pct,res_H1p5_err_pct,res_H1p25_err_pct], [deltas_H2,deltas_H1p5,deltas_H1p25],
                             '', ['$H/h=2$', '$H/h=1.5$', '$H/h=1.25$'],['$\delta$','$|\mathcal{R}_R-\mathcal{R}_S|/\mathcal{R}_S\%$'], loc='upper')

#------------------------------------------------------------------------------------------------------------
# dP_bfs = np.asarray([-40.12,-40.66,-41.55,-43.12,-51.65,-63.95])

# dP_cbfs = np.asarray([-40.01,-40.57,-41.48,-43.06,-51.58,-63.88])

# dP_pct_err = 100*(dP_cbfs-dP_bfs)/dP_bfs

# graphics.plot_2D(dP_pct_err, H[:-1], '', ['$H/h$','$|\Delta \mathcal{P}|\%$ error'])













