# -*- coding: utf-8 -*-
"""
Created on Tue Nov  5 08:26:55 2024

@author: sarah
"""
import numpy as np
import graphics

# BFS at Re=0, Q=2, U=0, L=4, varying H/h --> resistance error, reattachment lengths
# L=4
# h=1
# l=1

# H = np.asarray([2.75,2.5,2.25,2,1.5,1.25,1.125])
# dP = np.asarray([-40.12,-40.66,-41.55,-43.12,-51.65,-63.95,-75.91])

# flux_stokes = np.asarray([2,2,2,2,2,2,2])
# flux_reyn = np.asarray([2.92,2.84,2.74,2.61,2.28, 2.1, 2.04])

# res_stokes = dP/flux_stokes
# res_reyn = dP/flux_reyn
# res_err_pct = 100*np.abs(res_reyn-res_stokes)/np.abs(res_stokes)

# graphics.plot_2D_multi([res_reyn, res_stokes], H, 'Stokes vs Reynolds: Resistance', ['Reyn', 'Stokes'],['$H/h$','$R = dP/Q$'])
# graphics.plot_2D(res_err_pct, H, 'Stokes vs Reynolds: BFS resistance %-error', ['$H/h$','% error'])

#------------------------------------------------------------------------------------------------------------

# # dBFS at Re=0, Q=2, U=0, L=4, H=2, varying delta --> resistance error
H=2
h=1

deltas_H2=np.asarray([1.25, 1, 0.75, 0.5, 0.25, 1e-8])
deltas_H1p5=np.asarray([0.75, 0.5, 0.25, 0.125, 1e-8])

slopes_H2 = (2-h)/(2*deltas_H2)

slopes_H1p5 = (1.5-h)/(2*deltas_H1p5)

flux_stokes = 2

dP_H2=np.asarray([-44.52, -47.77, -50.85, -54.54, -58.87, -63.98])
flux_H2_reyn=np.asarray([2.08, 2.12, 2.15, 2.20, 2.28, 2.37])

dP_H1p5=np.asarray([-60, -62, -64.6, -66.6, -68.5])
flux_H1p5_reyn=np.asarray([2.04, 2.07,2.11,2.14,2.20])


res_H2_stokes = dP_H2/flux_stokes
res_H2_reyn = dP_H2/flux_H2_reyn
res_H1p5_stokes = dP_H1p5/flux_stokes
res_H1p5_reyn = dP_H1p5/flux_H1p5_reyn

res_H2_err_pct = 100*np.abs(res_H2_reyn-res_H2_stokes)/np.abs(res_H2_stokes)
res_H1p5_err_pct = 100*np.abs(res_H1p5_reyn-res_H1p5_stokes)/np.abs(res_H1p5_stokes)

graphics.plot_2D_inf(res_H2_err_pct, slopes_H2, '$\delta$-BFS Resistance error $H=2$', ['slope: $\Delta = (H-h)/(2\delta)$','% error'])
graphics.plot_2D_inf_multi([res_H2_stokes,res_H2_reyn], slopes_H2, 'Resistance $\delta$-BFS $H=2$', ['slope: $\Delta = (H-h)/(2\delta)$','R=dP/Q'],  ['stokes','reyn'])



graphics.plot_2D_inf(res_H1p5_err_pct, slopes_H1p5, '$\delta$-BFS Resistance error $H=1.5$', ['slope: $\Delta = (H-h)/(2\delta)$','% error'])
graphics.plot_2D_inf_multi([res_H1p5_stokes,res_H1p5_reyn], slopes_H1p5, 'Resistance $\delta$-BFS $H=1.5$', ['slope: $\Delta = (H-h)/(2\delta)$','R=dP/Q'],['stokes','reyn'])




#------------------------------------------------------------------------------------------------------------
# # no eddy BFS
# L=4
# h=1
# l=1

# H = np.asarray([2.75,2.5,2.25,2,1.5,1.25])
# dP_BFS = np.asarray([-40.12,-40.66,-41.55,-43.12,-51.65,-63.95])
# dP_noEddy = np.asarray([-40.01,-40.57,-41.48,-43.07,-51.58,-63.88])

# flux_stokes = np.asarray([2,2,2,2,2,2])
# flux_reyn_BFS = np.asarray([2.92,2.84,2.74,2.61,2.28,2.10])
# flux_reyn_noEddy = np.asarray([2.87,2.78,2.67,2.53,2.19,2.07])

# res_stokes_BFS = dP_BFS/flux_stokes
# res_reyn_BFS = dP_BFS/flux_reyn_BFS

# res_stokes_noEddy = dP_noEddy/flux_stokes
# res_reyn_noEddy = dP_noEddy/flux_reyn_noEddy

# res_err_pct_BFS = 100*np.abs(res_reyn_BFS-res_stokes_BFS)/np.abs(res_stokes_BFS)
# res_err_pct_noEddy = 100*np.abs(res_reyn_noEddy-res_stokes_noEddy)/np.abs(res_stokes_noEddy)

# res_err_pct_stokes=100*np.abs(res_stokes_BFS-res_stokes_noEddy)/np.abs(res_stokes_BFS)
# res_err_pct_reyn=100*np.abs(res_reyn_BFS-res_reyn_noEddy)/np.abs(res_reyn_BFS)
 
# x_r=np.asarray([0.47,0.43,0.40,0.35,0.25,0.15])
# y_r=np.asarray([0.51,0.48,0.45,0.41,0.3,0.17])

# # eddy_area=0.5*x_r*y_r
# # total_area= (h*l) + (H*(L-l))
# # eddy_pct_area=100*(1-(total_area-eddy_area)/total_area)

# graphics.plot_2D_multi([x_r,y_r], H, 'Re=0 flow separation points', ['$x_r$','$y_r$'], ['$H/h$','side lengths'])
# # graphics.plot_2D(eddy_pct_area, H, 'Re=0 BFS: eddy %-area', ['expansion ratio $H/h$','eddy %-area'])
 
# graphics.plot_2D_multi([res_reyn_BFS, res_stokes_BFS,res_reyn_noEddy, res_stokes_noEddy], H, 'Stokes vs Reynolds: Resistance', ['Reyn BFS', 'Stokes BFS', 'Reyn noEddy', 'Stokes noEddy'],['$H/h$','$R = dP/Q$'])

# graphics.plot_2D_multi([res_err_pct_BFS], H, 'Stokes vs Reynolds: resistance %-error',['BFS'], ['$H/h$','% error'])

# graphics.plot_2D_multi([res_err_pct_BFS,res_err_pct_noEddy], H, 'Stokes vs Reynolds: resistance %-error',['BFS','noEddy'], ['$H/h$','% error'])


# graphics.plot_2D_multi([res_err_pct_stokes], H, 'Effect of smoothing: resistance %-error',['stokes'], ['$H/h$','% error'])























