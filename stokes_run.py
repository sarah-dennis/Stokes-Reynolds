# -*- coding: utf-8 -*-
"""
Created on Tue Oct  1 11:56:12 2024

link to sd_run.sh

@author: sarah
"""

import stokes_control as control
import stokes_examples as examples
import graphics
import numpy as np

# example = examples.BFS_H2p5L4_Re0_Q2_U0
# example = examples.BFS_H2p25L4_Re0_Q2_U0
# example = examples.BFS_H2L4_Re0_Q2_U0
# example = examples.BFS_H1p5L4_Re0_Q2_U0
# example = examples.BFS_H1p25L4_Re0_Q2_U0
# example = examples.BFS_H1p125L4_Re0_Q2_U0

# H=[2.5,2.25,2,1.5,1.25,1.125]
# flux_err_BFS = [0.85,0.73,0.62,0.29,0.08,0.04]
# x_r=[0.43,0.39,0.36,0.25,0.156,0.085]
# y_r=[0.46,0.44,0.42,0.3,0.151,0.095]
# graphics.plot_2D(flux_err_BFS, H, 'Stokes vs Reynolds: Flux error', ['$H/h$','$Q_R-Q_S$'])
# graphics.plot_2D_multi([x_r,y_r], H, 'Re=0 attachment lengths', ['$x_r/h$','$y_r/h$'], ['$H/h$','length'])


#BLOWS UP (delta > 1.5, 1/3>slope)
# example = examples.BFS_H2L4_delta1p5 
# example = examples.BFS_H2L4_delta1
#BLOWS UP (1 > delta > 0.5, 1/2>slope>1)
# example = examples.BFS_H2L4_delta0p5 
# example = examples.BFS_H2L4_delta0p25
#BLOWS UP (0< delta <0.25, slope>2)
# example = examples.BFS_H2L4_delta0


#BLOWS UP (delta > 0.75, 1/3>slope)
# example = examples.BFS_H1p5L4_delta0p75
# example = examples.BFS_H1p5L4_delta0p5 
#BLOWS UP (0.5 > delta > 0.25, 1/2>slope>1)
# example = examples.BFS_H1p5L4_delta0p25
example = examples.BFS_H1p5L4_delta0p125 # running cluster
#BLOWS UP (0.125> delta >0, slope>2)
# example = examples.BFS_H1p5L4_delta0

# slopes=[1/3,1/2,1,2, np.inf]
# flux_err_dBFS_H2=[0.09, 0.12, 0.20,0.27, 0.36]
# flux_err_dBFS_H1p5=[0.02, 0.04, 0.09, 0.12, 0.19]
# graphics.plot_log_x_multi([flux_err_dBFS_H2,flux_err_dBFS_H1p5], slopes, 'Re=0 Stokes vs Reynolds Flux error', ['$H/h=2$','$H/h=1.5$'], ['slope: $\Delta = (H-h)/(2\delta)$','$Q_R-Q_S$'])



N = 200
#------------------------------------------------------------------------------
solver = control.Stokes_Solver(example) 

#-----------------------------------------------------------------------------

#------------------------------------------------------------------------------
# solver.new_run(N,1000)
# solver.new_run_many(N, 2, 3)
#------------------------------------------------------------------------------
# solver.load_scale(20,40)
solver.load_run(N,1)

#------------------------------------------------------------------------------
# solver.load_plot(N, zoom=False)
# ------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# BFS EXAMPLES  
#------------------------------------------------------------------------------

# TODO: BFS H=1.125, L= 4 Re=0, Q=2, U=0 --> N=200 re-run

# BFS H=2, L=4, Re=0, Q=2, U=0 --> done

# BFS H=1.5, L= 4 Re=0, Q=2, U=0 --> 640 needs more

# BFS H=1.25, L= 4 Re=0, Q=2, U=0 --> done

# solver.compare(20,[40,80,160,320],640)


#------------------------------------------------------------------------------
# BFS EDDIE REMOVAL EXAMPLES
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# Delta-BFS Examples
#------------------------------------------------------------------------------
