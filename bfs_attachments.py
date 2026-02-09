# -*- coding: utf-8 -*-
"""
Created on Fri Jan 30 13:25:51 2026

@author: sarah
"""

import numpy as np
import stokes_examples as examples
import graphics  
import readwrite as rw



def get_bfs_attachments(ex):
       
    u, v, psi, past_iters = rw.read_stokes(ex.filestr+".csv", ex.Nx * ex.Ny)
    
    y_xr = ex.yf   # xr on y=h
    x_yr = ex.xf/2 # yr on x=L/2
    
    i_yr =int(x_yr/ex.dx)-1 # 
    
    j_xr =int(y_xr/ex.dy)-1

    xrs = []
    yrs = []

    for i in range(ex.Nx-1):
        k_a = int(j_xr*ex.Nx + i)
        k_b = int(j_xr*ex.Nx + i+1 )   
        if np.sign(psi[k_a]-ex.flux)!= np.sign(psi[k_b]-ex.flux):
            
            xr = x_yr-(ex.xs[i]+ex.xs[i+1])/2
            if xr > 0:
                xrs.append(xr)
            
    for j in range(ex.Ny-1):
        k_a = int(j*ex.Nx + i_yr)
        k_b = int((j+1)*ex.Nx + i_yr)
        if np.sign(psi[k_a]-ex.flux)!= np.sign(psi[k_b]-ex.flux):
            yr = y_xr-(ex.ys[j] + ex.ys[j+1])/2
            if yr > 0:
                yrs.append(yr)
    if xrs==[]:
        xrs.append(0)
    if yrs == []:
        yrs.append(0)
        
    return xrs, yrs


U=0
Q=1
Re=0
#------------------------------------------------------------------------------

Example = examples.BFS
h_ins = [1.125, 1.25, 1.5, 2, 2.5, 2.75, 3] #N=80
# h_ins = [1.25, 2, 2.75] #N=160
h_out = 1
l_in = 8
l_out = 8

args_all = [[h_in, h_out, l_in, l_out] for h_in in h_ins]

num_tests = len(h_ins)

#------------------------------------------------------------------------------
# Example = examples.BFS_wedge

# h_in = 2
# h_out = 1
# l_in = 8
# l_out=8

# xyws = [[0.35,0.4],[0.2625, 0.3],[0.175,0.2]]

# args_all = [[h_in, h_out, l_in, l_out, xw, yw] for (xw, yw) in xyws]


# num_tests = len(xyws)


#------------------------------------------------------------------------------
xrs = np.zeros(num_tests)
yrs = np.zeros(num_tests)
N = 80

for k in range(num_tests):
    args = args_all[k]         
    ex = Example(args, U, Q, Re, N)
    xr, yr = get_bfs_attachments(ex)
    
    print(h_ins[k], xr, yr)
    #primary
    xrs[k] = xr[0]
    yrs[k] = yr[0]

print(xrs,yrs)

graphics.plot_2D_multi([xrs, yrs], h_ins, 'BFS Points of Flow Separation', ['$x_r$', '$y_r$'], ['$\mathcal{H}=H_{in}/H_{out}$', 'length'], loc='left')
    