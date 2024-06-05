# -*- coding: utf-8 -*-
"""
Created on Mon Jun  3 15:23:37 2024

@author: sarah
"""

import numpy as np

def l1pctError(f_u, f_v, g_u, g_v, Nx, Ny):
    max_err = 0.0
    norm_err = np.zeros((Ny, Nx))

    for j in range(Ny):
        for i in range(Nx):
                        
            fu = f_u[j,i]
            
            fv = f_v[j,i]
            
            gu = g_u[j,i]
            
            gv = g_v[j,i]
            
            err_ij = abs(fu - gu) + abs(fv - gv)
            
            norm_ij = abs(gu) + abs(gv)
            
            if norm_ij == 0:
                norm_err_ij = 0
            else:
                norm_err_ij = 100*err_ij/norm_ij
            norm_err[j,i] = norm_err_ij
            max_err = max(norm_err_ij, max_err)
    return norm_err, max_err

def l2pctError(f_u, f_v, g_u, g_v, Nx, Ny):
    max_err = 0.0
    
    norm_err = np.zeros((Ny, Nx))
    for j in range(Ny):
        for i in range(Nx):
            fu = f_u[j,i]
            
            fv = f_v[j,i]
            
            gu = g_u[j,i]
            
            gv = g_v[j,i]
           
            err_ij = np.sqrt((fu - gu)**2 + (fv - gv)**2)

            norm_ij = np.sqrt(gu**2 + gv**2)
            
            if norm_ij == 0:
                norm_err_ij = 0
            else:
                norm_err_ij = 100* err_ij/norm_ij
          
            norm_err[j,i] = norm_err_ij
            
            max_err = max(norm_err_ij, max_err)
    return norm_err, max_err

def lInfpctError(f_u, f_v, g_u, g_v, Nx, Ny):

    max_err = 0.0
    norm_err = np.zeros((Ny, Nx))

    for j in range(Ny):
        for i in range(Nx):
            fu = f_u[j,i]
            
            fv = f_v[j,i]
            
            gu = g_u[j,i]
            
            gv = g_v[j,i]
           
            err_ij = max(abs(fu - gu), abs(fv - gv))
            
            norm_ij = max(abs(gu), abs(gv))
            
            if norm_ij == 0:
                norm_err_ij = 0
            else:
                norm_err_ij = 100* err_ij/norm_ij
            norm_err[j,i] = norm_err_ij
            max_err = max(norm_err_ij, max_err)
            
    return norm_err, max_err

