# -*- coding: utf-8 -*-
"""
Created on Fri Jan 30 13:25:51 2026

@author: sarah
"""

import numpy as np
import stokes_examples as examples
import graphics  
import readwrite as rw



def get_tri_attachments(ex, H, L):
       
    u, v, psi, past_iters = rw.read_stokes(ex.filestr+".csv", ex.Nx * ex.Ny)
    xrs=[]
    for i in range(ex.Nx):
        if i < ex.Nx/2 - 1:
            j_ha = H/(L/2) * i
            j_hb = H/(L/2) * (i+1)
            k_a = int(j_ha*ex.Nx + i)
            k_b = int(j_hb*ex.Nx + i+1)   
            
          
            if np.sign(psi[k_a]-ex.flux)!= np.sign(psi[k_b]-ex.flux):
                
                xr = ex.xs[i]
                xrs.append(xr)


    if xrs==[]:
        xrs.append(0)

        
    return xrs


U=1
Q=0
Re=0
#------------------------------------------------------------------------------


hs = [4, 2, 1]
L = 2
N=160

#------------------------------------------------------------------------------

num_h = len(hs)


for i in range(num_h):
    H= hs[i]

    Example = examples.TriCavity
    args = [H, L]           
    ex = Example(args, U, Q, Re, N)
    
    xr = get_tri_attachments(ex, H, L)
    yr = [2*H/L * x for x in xr]
    print(H, xr, yr)


    