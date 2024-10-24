# -*- coding: utf-8 -*-
"""
Created on Fri May 24 08:36:17 2024

@author: sarah
"""
import numpy as np
import stokes_readwrite as rw
#-------------------------------------------------------------------------------

def compare_Ns(ex, N_min, Ns, N_max):
    ex_min = ex(N_min)

    ex_max = ex(N_max)
    u_max, v_max, psi_max, past_iters = rw.read_solution(ex_max.filestr+".csv", ex_max.Nx * ex_max.Ny)
    psi_max = psi_max.reshape((ex_max.Ny,ex_max.Nx))
    mult_max = int(N_max/N_min)
    err = np.zeros((len(Ns)+1,ex_min.Ny,ex_min.Nx))
    err_inf = np.zeros(len(Ns)+1)
    err_l1 = np.zeros(len(Ns)+1)
    err_l2 = np.zeros(len(Ns)+1)
    
    for n in range(len(Ns)+1):
        if n == 0:
            ex_n = ex(N_min)
            u, v, psi_n, past_iters = rw.read_solution(ex_n.filestr+".csv", ex_n.Nx * ex_n.Ny)
            psi_n=psi_n.reshape((ex_min.Ny, ex_min.Nx))
            
            
            mult = 1
        else:
            ex_n = ex(Ns[n-1])
            u, v, psi_n, past_iters = rw.read_solution(ex_n.filestr+".csv", ex_n.Nx * ex_n.Ny)
            psi_n=psi_n.reshape((ex_n.Ny, ex_n.Nx))
            mult = int(Ns[n-1]/N_min)

        for k_min in range(ex_min.Ny*ex_min.Nx):
        # all indices (i,j) on grid N_min
            i = k_min % ex_min.Nx
            j = (k_min // ex_min.Nx)
            
            # -> indices on grid N || N_min
            i_n = mult * i
            j_n = mult * j

            # -> indices on grid N_max || N_min
            i_max = mult_max * i
            j_max = mult_max * j
            
            err[n,j,i] = abs(psi_max[j_max,i_max] - psi_n[j_n,i_n])

            err_l1[n] += err[n,j,i]
            err_l2[n] += err[n,j,i]**2
        # graphics.plot_contour_mesh(err[n], ex_min.xs, ex_min.ys, "err", ['x', 'y', 'err'], True)   
        err_inf[n] = np.max(err[n])
        err_l2[n] = np.sqrt(err_l2[n])

    l1_rate = convg_rate(err_l1)
    l2_rate = convg_rate(err_l2)
    inf_rate = convg_rate(err_inf)
    cnvg_rates = np.stack([l1_rate, l2_rate, inf_rate], axis=0)

    return err_l1, err_l2, err_inf, cnvg_rates, ex_min
        
def convg_rate(errs):
    n = len(errs)
    rates = np.zeros(n-1)
    for k in range(n-1):
        rates[k]=errs[k+1]/errs[k]
    
    return rates







#------------------------------------------------------------------------------
# For BFS examples
#------------------------------------------------------------------------------
def get_attatchments(ex, N):

    step = ex(N)
    u, v, psi, past_iters = rw.read_solution(step.filestr+".csv", step.Nx * step.Ny)
    
    psi_2D = psi.reshape((step.Ny,step.Nx))
    psi_xs_yf = psi_2D[step.jf_out] #reattatchment wall 
    psi_ys_xstep = psi_2D[:,step.i_step] # detatchment wall
    
    xs_saddle = []
    sign_ref = 0
    for i in range(1, step.Nx):
        sign_new = np.sign(psi_xs_yf[i])
        if sign_new != sign_ref:
            sign_ref = sign_new
            xs_saddle.append(step.xs[i])

    
    ys_saddle = []
    sign_ref = 0
    for j in range(1, step.Ny):
        sign_new = np.sign(psi_ys_xstep[j])
        if sign_new != sign_ref:
            sign_ref = sign_new
            ys_saddle.append(step.ys[j])

        
    return xs_saddle, ys_saddle
    

#------------------------------------------------------------------------------


