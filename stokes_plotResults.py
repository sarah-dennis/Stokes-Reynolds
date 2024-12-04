# -*- coding: utf-8 -*-
"""
Created on Tue Nov  5 08:26:55 2024

@author: sarah
"""
import numpy as np
import graphics

# BFS at Re=0, Q=2, U=0, L=4, varying H/h --> resistance error, reattachment lengths
# h=1
# L=4
# l=1
# H = np.asarray([2.75,2.5,2.25,2,1.5,1.25,1.125])

# dP_stokes = np.asarray([-40.12,-40.66,-41.55,-43.12,-51.65,-63.95,-75.91])
# dP_reyn = np.asarray([-27.50,-28.55,-30.25,-33.05,-45.25,-61.00,-74.50])
# flux = np.asarray([2,2,2,2,2,2,2])
# # flux_reyn = np.asarray([2.92,2.84,2.74,2.61,2.28, 2.1, 2.04])

# res_stokes = dP_stokes/flux
# res_reyn = dP_reyn/flux
# res_err_pct = 100*np.abs(res_reyn-res_stokes)/np.abs(res_stokes)

# print('Resistance %err:', res_err_pct)

# graphics.plot_2D_multi([res_reyn, res_stokes], H, 'Stokes vs Reynolds: BFS Resistance', ['Reyn', 'Stokes'],['$H/h$','$R = \Delta P / Q$'])
# graphics.plot_2D(res_err_pct, H, 'Stokes vs Reynolds: BFS resistance %-error', ['$H/h$','Resistance $R$ % error'], color='darkviolet')

# x_r=np.asarray([0.47,0.43,0.40,0.35,0.25,0.15,0.09])
# y_r=np.asarray([0.51,0.48,0.45,0.41,0.3,0.17,0.1])

# eddy_area=0.5*x_r*y_r
# total_area= (h*l) + (H*(L-l))
# eddy_pct_area=100*(1-(total_area-eddy_area)/total_area)

# graphics.plot_2D_multi([x_r,y_r], H, 'Re=0 BFS Flow Separation Points', ['$x_r$','$y_r$'], ['$H$','lengths $|x_r|$, $|y_r|$'])
# graphics.plot_2D(eddy_pct_area, H, 'Re=0 BFS: eddy %-area', ['expansion ratio $H/h$','eddy %-area'])
 



#------------------------------------------------------------------------------------------------------------

# # dBFS at Re=0, Q=2, U=0, L=4, H=2, varying delta --> resistance error

deltas_H2=np.asarray([1, 0.75, 0.5, 0.25, 0.125, 0])
deltas_H1p5=np.asarray([0.75, 0.5, 0.25, 0.125, 0])
deltas_H1p25=np.asarray([0.25, 0.125, 0.05, 0])

flux_stokes = 2

dP_H2=np.asarray([-47.8, -50.85, -54.54, -58.87,-61.8, -63.98])
flux_H2_reyn=np.asarray([2.12, 2.15, 2.20, 2.28, 2.34, 2.37])

dP_H1p5=np.asarray([-59.89, -61.98, -64.58, -66.64, -68.46])
flux_H1p5_reyn=np.asarray([2.04, 2.07, 2.11, 2.16, 2.20])


dP_H1p25=np.asarray([-73.2, -74.2, -75, -75.6])
flux_H1p25_reyn=np.asarray([2.03,2.05,2.07,2.08])

res_H2_stokes = dP_H2/flux_stokes
res_H2_reyn = dP_H2/flux_H2_reyn
res_H1p5_stokes = dP_H1p5/flux_stokes
res_H1p5_reyn = dP_H1p5/flux_H1p5_reyn
res_H1p25_stokes = dP_H1p25/flux_stokes
res_H1p25_reyn = dP_H1p25/flux_H1p25_reyn

graphics.plot_2D_multi([res_H2_stokes,res_H2_reyn], deltas_H2, 'Resistance $\delta$-BFS \n $H=2$, $L=4, Re=0, Q=2, U=0$', ['stokes','reyn'], ['$\delta$','R=dP/Q'])
graphics.plot_2D_multi([res_H1p5_stokes,res_H1p5_reyn], deltas_H1p5, 'Resistance $\delta$-BFS \n $H=1.5$, $L=4, Re=0, Q=2, U=0$', ['stokes','reyn'], ['$\delta$','R=dP/Q'])
graphics.plot_2D_multi([res_H1p25_stokes,res_H1p25_reyn], deltas_H1p25, 'Resistance $\delta$-BFS \n $H=1.25$, $L=4, Re=0, Q=2, U=0$', ['stokes','reyn'], ['$\delta$','R=dP/Q'])

graphics.plot_2D_multi_multi([res_H2_reyn,res_H2_stokes,res_H1p5_reyn,res_H1p5_stokes,res_H1p25_reyn,res_H1p25_stokes],[deltas_H2,deltas_H2,deltas_H1p5,deltas_H1p5,deltas_H1p25,deltas_H1p25],
                              'Resistance $\delta$-BFS \n $L=4, Re=0, Q=2, U=0$',  ['$H=2$ reyn','$H=2$ stokes','$H=1.5$ reyn','$H=1.5$ stokes','$H=1.25$ reyn','$H=1.25$ stokes'], ['$\delta$','R=dP/Q'])

res_H2_err_pct = 100*np.abs(res_H2_reyn-res_H2_stokes)/np.abs(res_H2_stokes)
res_H1p5_err_pct = 100*np.abs(res_H1p5_reyn-res_H1p5_stokes)/np.abs(res_H1p5_stokes)
res_H1p25_err_pct = 100*np.abs(res_H1p25_reyn-res_H1p25_stokes)/np.abs(res_H1p25_stokes)

graphics.plot_2D_multi_multi([res_H2_err_pct,res_H1p5_err_pct,res_H1p25_err_pct], [deltas_H2,deltas_H1p5,deltas_H1p25], 'Reynolds to Stokes Resistance error \n $\delta$-BFS, $L=4, Re=0, Q=2, U=0$', ['$H=2$', '$H=1.5$', '$H=1.25$'],['$\delta$','% error: $(R_{r} - R_{s})/R_{s}$'])



















