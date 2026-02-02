# -*- coding: utf-8 -*-
"""
Created on Fri Jan 30 13:25:51 2026

@author: sarah
"""

import numpy as np
import stokes_control as control
import stokes_examples as examples
import graphics  

U=0
Q=1
Re=0
#------------------------------------------------------------------------------

# Example = examples.BFS

# h_ins = [1.125, 1.25, 1.5, 2, 2.5, 2.75, 3]
# h_out = 1
# l_in = 8
# l_out = 8

# args_all = [[h_in, h_out, l_in, l_out] for h_in in h_ins]

# num_tests = len(h_ins)

#------------------------------------------------------------------------------
Example = examples.BFS_wedge

h_in = 2
h_out = 1
l_in = 8
l_out=8

xyws = [[0.35,0.4],[0.2625, 0.3],[0.175,0.2]]

args_all = [[h_in, h_out, l_in, l_out, xw, yw] for (xw, yw) in xyws]


num_tests = len(xyws)


#------------------------------------------------------------------------------
xrs = np.zeros(num_tests)
yrs = np.zeros(num_tests)
N = 160

for k in range(num_tests):
    args = args_all[k]
    solver = control.Stokes_Solver(Example, args, U, Q, Re, max_iters=500000)                

    xr, yr = solver.get_bfs_attachments(N)
    
    
    #primary
    xrs[k] = xr[0]
    yrs[k] = yr[0]
print(xrs,yrs)


# graphics.plot_2D_multi([xrs, yrs], h_ins, 'BFS Flow Stagnation Points', ['$x_r$', '$y_y$'], ['$\mathcal{H}=H_{in}/H_{out}$', 'length'], loc='left')
    