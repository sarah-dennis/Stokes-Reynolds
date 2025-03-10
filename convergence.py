# -*- coding: utf-8 -*-
"""
Created on Fri May 24 08:36:17 2024

@author: sarah
"""
import numpy as np
import readwrite as rw
import stokes_pressure as pressure
#-------------------------------------------------------------------------------
def stokes_cnvg_self(Ex, N_min, Ns, N_max, p_err=True):
    ex_min = Ex(N_min)
    ex_max = Ex(N_max)
    
    # Load max grid for 'true' values
    u_max, v_max, psi_max, past_iters = rw.read_stokes(ex_max.filestr+".csv", ex_max.Nx * ex_max.Ny)
    psi_max = psi_max.reshape((ex_max.Ny,ex_max.Nx))

    if p_err:
        p_max = pressure.pressure(ex_max, u_max, v_max)
        dp_max, res_max = pressure.resistance(ex_max, p_max) 
        p_max = p_max.reshape((ex_max.Ny,ex_max.Nx))
        
    mult_max = int(N_max/N_min)
    
    err = np.zeros((len(Ns)+1,ex_min.Ny,ex_min.Nx))
    err_inf = np.zeros(len(Ns)+1)
    err_l1 = np.zeros(len(Ns)+1)
    err_l2 = np.zeros(len(Ns)+1)
    

    if p_err:
        err_p = np.zeros((len(Ns)+1,ex_min.Ny,ex_min.Nx))
        err_p_inf = np.zeros(len(Ns)+1)
        err_p_l1 = np.zeros(len(Ns)+1)
        err_p_l2 = np.zeros(len(Ns)+1)
        err_dp = np.zeros(len(Ns)+1)
    
    for n in range(len(Ns)+1):
        if n == 0:
            ex_n = Ex(N_min)
            mult = 1
        else:
            ex_n = Ex(Ns[n-1])
            mult = int(Ns[n-1]/N_min)
            
            
        u_n, v_n, psi_n, past_iters = rw.read_stokes(ex_n.filestr+".csv", ex_n.Nx * ex_n.Ny)
        psi_n=psi_n.reshape((ex_n.Ny, ex_n.Nx))

        if p_err:
            p_n = pressure.pressure(ex_n, u_n, v_n)
            dp_n, res_n = pressure.resistance(ex_n, p_n)
            p_n=p_n.reshape((ex_n.Ny,ex_n.Nx))
            
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
            
            if p_err:
                err_p[n,j,i] = abs(p_max[j_max,i_max] - p_n[j_n,i_n])
                err_p_l1[n]+= err_p[n,j,i] 
                err_p_l2[n]+= (err_p[n,j,i]**2)
                
        size_n= ex_n.Nx*ex_n.Ny   
        
        err_inf[n] = np.max(err[n])
        err_l1[n] /= size_n
        err_l2[n] = np.sqrt(err_l2[n] / size_n)
        
        if p_err:
            err_p_inf[n] = np.max(err_p[n])
            err_p_l2[n] = np.sqrt(err_p_l2[n])
            err_dp[n]=np.abs(dp_max-dp_n)
            
    l1_rate = convg_rate(err_l1)
    l2_rate = convg_rate(err_l2)
    inf_rate = convg_rate(err_inf)
    cnvg_rates = np.stack([l1_rate, l2_rate, inf_rate], axis=0)
    
    if p_err:
        p_l1_rate = convg_rate(err_p_l1)
        p_l2_rate = convg_rate(err_p_l2)
        p_inf_rate = convg_rate(err_p_inf)
        dp_rate = convg_rate(err_dp)
        p_cnvg_rates = np.stack([p_l1_rate, p_l2_rate, p_inf_rate, dp_rate], axis=0)
        
    print("stream cnvg rates")
    print("l1: " + np.array2string(cnvg_rates[0], precision=2))
    print("l2: " + np.array2string(cnvg_rates[1], precision=2))
    print("linf" + np.array2string(cnvg_rates[2], precision=2))
    
    if p_err:
        print("pressure cnvg rates")
        print("l1: " + np.array2string(p_cnvg_rates[0], precision=2))
        print("l2: " + np.array2string(p_cnvg_rates[1], precision=2))
        print("linf" + np.array2string(p_cnvg_rates[2], precision=2))
        print("dp: " + np.array2string(p_cnvg_rates[3], precision=2))

    return err_l1, err_l2, err_inf, cnvg_rates, ex_min

#------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
def reyn_cnvg_self(solver, N_min, Ns, N_max):
    # Load max grid for 'true' values
    p_max, vel_max = solver.fd_adj_solve(N_max, plot=True)
    p_min, vel_min = solver.fd_adj_solve(N_min, plot=False)
    
    ps_max = np.nan_to_num(p_max.ps_2D)
    ps_min = np.nan_to_num(p_min.ps_2D)

    Nx_min = ps_min.shape[1]
    Ny_min = ps_min.shape[0]


    us_max = np.nan_to_num(vel_max.vx)
    vs_max = np.nan_to_num(vel_max.vy)
        

    mult_max = int(N_max/N_min)
    

    err_p = np.zeros((len(Ns)+1,Ny_min,Nx_min))
    err_p_inf = np.zeros(len(Ns)+1)
    err_p_l1 = np.zeros(len(Ns)+1)
    err_p_l2 = np.zeros(len(Ns)+1)
    
    err_vx = np.zeros((len(Ns)+1,Ny_min,Nx_min))
    err_vy = np.zeros((len(Ns)+1,Ny_min,Nx_min))
    err_vel_inf = np.zeros(len(Ns)+1)
    err_vel_l1 = np.zeros(len(Ns)+1)
    err_vel_l2 = np.zeros(len(Ns)+1)
    
    for n in range(len(Ns)+1):
        if n == 0:
            p_n, vel_n = solver.fd_adj_solve(N_min, plot=False)
            mult = 1
        else:
            p_n, vel_n = solver.fd_adj_solve(Ns[n-1], plot=False)
        
            mult = int(Ns[n-1]/N_min)
        
        ps_n = np.nan_to_num(p_n.ps_2D)
        us_n = np.nan_to_num(vel_n.vx)
        vs_n = np.nan_to_num(vel_n.vy)
            
        for k_min in range(Ny_min*Nx_min):
        # all indices (i,j) on grid N_min
            i = k_min % Nx_min
            j = (k_min // Nx_min)
            
            # -> indices on grid N || N_min
            i_n = mult * i
            j_n = mult * j

            # -> indices on grid N_max || N_min
            i_max = mult_max * i
            j_max = mult_max * j
            
            err_p[n,j,i] = abs(ps_max[j_max,i_max] - ps_n[j_n,i_n])
            err_p_l1[n] += err_p[n,j,i]
            err_p_l2[n] += err_p[n,j,i]**2

            err_vx[n,j,i] = abs(us_max[j_max,i_max] - us_n[j_n,i_n])
            err_vy[n,j,i] = abs(vs_max[j_max,i_max] - vs_n[j_n,i_n])
            err_vel_l1[n]+= err_vx[n,j,i] + err_vy[n,j,i]
            err_vel_l2[n]+= (err_vx[n,j,i]**2 + err_vy[n,j,i]**2)
                
        size_n=  ps_n.shape[1]*ps_n.shape[0]
        
        err_p_inf[n] = np.max(err_p[n])
        err_p_l1[n] /= size_n
        err_p_l2[n] = np.sqrt(err_p_l2[n] / size_n)
        
    
        err_vel_inf[n] = np.max([np.max(err_vx[n]),np.max(err_vy[n])])
        err_vel_l1[n] /= size_n
        err_vel_l2[n] = np.sqrt(err_vel_l2[n]/size_n)

            
    p_l1_rate = convg_rate(err_p_l1)
    p_l2_rate = convg_rate(err_p_l2)
    p_inf_rate = convg_rate(err_p_inf)
    p_cnvg_rates = np.stack([p_l1_rate, p_l2_rate, p_inf_rate], axis=0)
    

    vel_l1_rate = convg_rate(err_vel_l1)
    vel_l2_rate = convg_rate(err_vel_l2)
    vel_inf_rate = convg_rate(err_vel_inf)
    vel_cnvg_rates = np.stack([vel_l1_rate, vel_l2_rate, vel_inf_rate], axis=0)
    
    print("pressure cnvg rates")
    print("l1: " + np.array2string(p_cnvg_rates[0], precision=2))
    print("l2: " + np.array2string(p_cnvg_rates[1], precision=2))
    print("linf" + np.array2string(p_cnvg_rates[2], precision=2))
    

    print("velocity cnvg rates")
    print("l1: " + np.array2string(vel_cnvg_rates[0], precision=2))
    print("l2: " + np.array2string(vel_cnvg_rates[1], precision=2))
    print("linf" + np.array2string(vel_cnvg_rates[2], precision=2))

    p_info = [err_p_l1, err_p_l2, err_p_inf, p_cnvg_rates]
    vel_info = [err_vel_l1, err_vel_l2, err_vel_inf, vel_cnvg_rates]
    return p_info, vel_info

#------------------------------------------------------------------------------



def reyn_cnvg_pwl_fd(solver, Ns):
    many=len(Ns)
    inf_errs = np.zeros(many)
    l1_errs = np.zeros(many)
    l2_errs = np.zeros(many)
    q_errs = np.zeros(many)


    for k in range (many):
        N=Ns[k]
        pwl_p, pwl_vel = solver.pwl_solve(N, plot=False)
        fd_p, fd_vel = solver.fd_solve(N, plot=False)
        size_k = len(pwl_p.ps_1D)
        
        p_err = np.abs(pwl_p.ps_1D - fd_p.ps_1D)
        inf_errs[k] = np.max(p_err)
        l1_errs[k] = np.sum(p_err)/size_k
        l2_errs[k] = np.sqrt(np.sum(p_err**2)/size_k)
        q_errs[k] = np.abs(pwl_vel.flux-fd_vel.flux)
        
    p_l1_rate = convg_rate(l1_errs)
    p_l2_rate = convg_rate(l2_errs)
    p_inf_rate = convg_rate(inf_errs)
    q_err_rate = convg_rate(q_errs)
    cnvg_rates = np.stack([p_l1_rate, p_l2_rate, p_inf_rate, q_err_rate], axis=0)
    
    
        
    print("cnvg rates")
    print("p-l1: " + np.array2string(p_l1_rate))
    print("p-l2: " + np.array2string(p_l2_rate))
    print("p-linf" + np.array2string(p_inf_rate))
    print("q-linf" + np.array2string(q_err_rate))
    
    return l1_errs, l2_errs, inf_errs, cnvg_rates
    
def convg_rate(errs):
    n = len(errs)
    rates = np.zeros(n-1)
    
    for k in range(n-1):
        rates[k]=errs[k+1]/errs[k]
    
    return rates

#------------------------------------------------------------------------------


