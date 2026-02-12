# -*- coding: utf-8 -*-
"""
Created on Fri Jan 30 13:25:51 2026

@author: sarah
"""

import numpy as np
import stokes_control as control
import stokes_examples as examples
import graphics  
import readwrite as rw



def get_bfsPWL_attachments(ex, delta):
       
    u, v, psi, past_iters = rw.read_stokes(ex.filestr+".csv", ex.Nx * ex.Ny)
    
    y_xr = ex.yf   # xr on y=h  
    j_xr =int(y_xr/ex.dy)-1
    
    x_corner = (ex.xf - delta)/2

    xrs = []

    for i in range(ex.Nx-1):
        k_a = int(j_xr*ex.Nx + i)
        k_b = int(j_xr*ex.Nx + i+1)   
        if np.sign(psi[k_a]-ex.flux)!= np.sign(psi[k_b]-ex.flux):
            
            xr = x_corner-(ex.xs[i]+ex.xs[i+1])/2
            if xr > 0:
                xrs.append(xr)

    if xrs==[]:
        xrs.append(0)

        
    return xrs


U=0
Q=1
Re=0
#------------------------------------------------------------------------------


h_ins = [ 2.75,2, 1.25]
deltas = [1, 1/2, 1/4, 1/8, 0]
h_out = 1
L=16

N=80

#------------------------------------------------------------------------------

num_h = len(h_ins)
num_d = len(deltas)
num_tests = num_d*num_h

xrs = np.zeros((num_h, num_d))
for h in range(num_h):
    for d in range(num_d):
        h_in = h_ins[h]
        delta = deltas[d]
        
        if delta != 0:
            Example = examples.BFS_pwl
            args = [h_in, h_out, L, delta]           
            ex = Example(args, U, Q, Re, N)
            xr = get_bfsPWL_attachments(ex, delta)
            xrs[h,d] = xr[0]
        else:
            Example = examples.BFS
            args = [h_in, h_out, L//2, L//2]           
            ex = Example(args, U, Q, Re, N)
            
            xr = get_bfsPWL_attachments(ex, delta)
            xrs[h,d] = xr[0]
        
print(xrs)

fun_labels = ['$\mathcal{H}=%.2f$'%h_in for h_in in h_ins]
graphics.plot_2D_multi(xrs, deltas, 'Wedge Slider Points of Flow Separation $x_r$', fun_labels, ['$\delta$', 'length'], loc='upper')

    