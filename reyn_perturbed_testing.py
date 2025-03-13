# -*- coding: utf-8 -*-
"""
Created on Wed Mar  5 12:53:02 2025

@author: sarah
"""


import reyn_control as control
import reyn_examples as examples
import convergence

import numpy as np
import graphics

#------------------------------------------------------------------------------
# Convergence to analytic solution
#------------------------------------------------------------------------------
Example = examples.LambdaBump # 
H=1 #=h0
l=1 #=L/2
delta = H/l
lam = 0.8
args=[lam, H, l]

U=0
dP=1

n_Ns = 5
Ns = np.zeros(n_Ns)
errs = np.zeros(n_Ns)
errs_p0 = np.zeros(n_Ns)
errs_p2 = np.zeros(n_Ns)
dN = 2
N = 7

# analytic solution dP for LambdaBump
ardP = -3*(3*(lam**2)-8*lam+8)/(1-lam)**(5/2)
padP = -12*((3.14159*(H/l)*lam)**2)/(5*(1-lam)**(3/2))
adP = (ardP + padP) 
 
for k in range(n_Ns):
    Ns[k] = N 
    
    solver = control.Reynolds_Solver(Example, U, dP, args)
    rdP, p2dP, p2p, p2v = solver.fd_pert_solve(N, write=False, plot=False, get_dPs=True)

    errs[k] = abs(p2dP -adP)/abs(adP)
    errs_p0[k] = abs(rdP -ardP)/abs(ardP)
    errs_p2[k] = abs((p2dP-ardP) -padP)/abs(padP)
    
    N *= dN
    
linthresh = min(errs)*1e-1
graphics.plot_log_multi([errs, errs_p0, errs_p2], Ns, f'ELT-2 convergence $(\lambda={lam}, \delta={delta})$',[ '$\Delta P$' , '$\Delta P_0$', '$\Delta P_2$'], ['N', 'Rel. Error '], linthresh, O1=1, O2=1e1)


p_rate = convergence.convg_rate(errs)
p0_rate = convergence.convg_rate(errs_p0)
p2_rate = convergence.convg_rate(errs_p2)
cnvg_rates = np.stack([p_rate, p0_rate, p2_rate], axis=0)


    
print("cnvg rates")
print("p: " + np.array2string(p_rate))
print("p0: " + np.array2string(p0_rate))
print("p2" + np.array2string(p2_rate))

# #------------------------------------------------------------------------------
# # Comparison of dP with Reynolds (unadjusted) model
# #------------------------------------------------------------------------------
# Example = examples.LambdaBump # can change to any example
# H=1 #=h0
# l=1 #=L/2
# delta = H/l
# lam = 0.8
# args=[lam, H, l]

# U=0
# dP=1

# N=100
# n_lams =5
# dlam= 1/(n_lams+2)
# pdPs = np.zeros(n_lams)
# rdPs = np.zeros(n_lams)
# adPs = np.zeros(n_lams)
# lams = np.zeros(n_lams)
# for k in range(n_lams):
#     lam= dlam * (k+1)

#     args=[lam, H, l]
#     solver = control.Reynolds_Solver(Example, U, dP, args)
#     rdP, pdP, pp, pv = solver.fd_pert_solve(N, write=False, plot=False, get_dPs=True)

#     rp, rv = solver.fd_solve(N, write=False, plot=False)
#     lams[k] = lam
#     rdPs[k] = -rdP
#     pdPs[k] = -pdP
    
#     ardP = 3*(3*lam**2-8*lam+8)/(1-lam)**(5/2)
#     padP = 12*(np.pi*(H/l)*lam)**2/(5*(1-lam)**(3/2))
#     adPs[k] = ardP + padP 

# graphics.plot_2D_multi([rdPs,adPs,pdPs], lams, f'CLT vs ELT-2 $(\delta={delta})$', ['CLT', 'analytic ELT-2','ELT-2'], ['$\lambda$','$\Delta P$'], loc='left')
# graphics.plot_2D(100*abs(adPs-pdPs)/abs(pdPs), lams, '%-err ELT-analytic vs ELT-2 (arch)',  ['$\lambda$','err $\Delta P$'])
