# -*- coding: utf-8 -*-
"""
Created on Tue Nov  5 08:26:55 2024

@author: sarah
"""
import numpy as np
import graphics

# BFS at Re=0, Q=2, U=0, L=4, varying H/h --> resistance error, reattachment lengths
L=4
h=1
l=1

H = np.asarray([2.75,2.5,2.25,2,1.5,1.25,1.125])
dP = np.asarray([-52.81,-52.44,-51.62,-50.11,-42.21,-31.82,-24.67])

flux_stokes = np.asarray([2,2,2,2,2,2,2])
flux_reyn = np.asarray([3.85,3.67,3.40, 3.04,1.86, 1.05, 0.66])

res_stokes = dP/flux_stokes
res_reyn = dP/flux_reyn
res_err_pct = 100*np.abs(res_reyn-res_stokes)/np.abs(res_stokes)

x_r=np.asarray([0.43,0.43,0.39,0.36,0.25,0.151,0.085])
y_r=np.asarray([0.5,0.46,0.44,0.42,0.3,0.153,0.095])
eddy_area=0.5*x_r*y_r
total_area= (h*l) + (H*(L-l))
eddy_pct_area=100*(1-(total_area-eddy_area)/total_area)

graphics.plot_2D(res_err_pct, H, 'Stokes vs Reynolds: BFS resistance %-error', ['$H/h$','% error'])

graphics.plot_2D_multi([res_reyn, res_stokes], H, 'Stokes vs Reynolds: Resistance', ['Reyn', 'Stokes'],['$H/h$','$R = dP/Q$'])

graphics.plot_2D_multi([x_r,y_r], H, 'Re=0 flow separation points', ['$x_r$','$y_r$'], ['$H/h$','side lengths'])

graphics.plot_2D(eddy_pct_area, H, 'Re=0 BFS: eddy %-area', ['expansion ratio $H/h$','eddy %-area'])

#------------------------------------------------------------------------------------------------------------


# dBFS at Re=0, Q=2, U=0, L=4, H=2, varying delta --> resistance error

# slopes=[1/3, 1/2, 1, 2, 1000] 
# flux_stokes = np.asarray([2,2,2,2,2])

# dP_H2=np.asarray([-42.35,-47.69,-54.41,-58.77,-63.83])
# flux_H2_reyn=np.asarray([2.09, 2.12, 2.20, 2.27, 2.36])

# dP_H1p5=np.asarray([-59.47,-61.23,-63.89,-64.5,-68.23])
# flux_H1p5_reyn=np.asarray([2.02, 2.04, 2.09, 2.12, 2.19])

# res_H2_stokes = dP_H2/flux_stokes
# res_H2_reyn = dP_H2/flux_H2_reyn
# res_H1p5_stokes = dP_H1p5/flux_stokes
# res_H1p5_reyn = dP_H1p5/flux_H1p5_reyn
# res_H2_err_pct = 100*np.abs(res_H2_reyn-res_H2_stokes)/np.abs(res_H2_stokes)
# res_H1p5_err_pct = 100*np.abs(res_H1p5_reyn-res_H1p5_stokes)/np.abs(res_H1p5_stokes)

# graphics.plot_2D_inf_multi([res_H2_stokes,res_H1p5_stokes], slopes, 'Stokes resistance ', ['slope: $\Delta = (H-h)/(2\delta)$','R=dP/Q'], ['H=2','H=1.5'])

# graphics.plot_2D_inf_multi([res_H2_reyn,res_H1p5_reyn], slopes, 'Reynolds resistance ', ['slope: $\Delta = (H-h)/(2\delta)$','R=dP/Q'], ['H=2','H=1.5'])

# graphics.plot_2D_inf_multi([res_H2_err_pct,res_H1p5_err_pct], slopes, 'Stokes vs Reynolds resistance %-error', ['slope: $\Delta = (H-h)/(2\delta)$','% error'], ['H=2','H=1.5'])
