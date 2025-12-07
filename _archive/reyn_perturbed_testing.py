# -*- coding: utf-8 -*-
"""
Created on Wed Mar  5 12:53:02 2025

@author: sarah
"""


import reyn_control as control
import reyn_examples as examples
import convergence
import reyn_boundary as rbc

import numpy as np
import graphics

#------------------------------------------------------------------------------
# Convergence to analytic solution
#------------------------------------------------------------------------------
Example = examples.LambdaBump # 
H=1 #=h0
l=1 #=L/2
delta = H/l
lam = 0.9 # < 1
h0=1
args=[lam, H, l,h0]

# lam=1/2
# H=2
# l=4
# h0 = 3/2
# args=[lam, H, l, h0]

U=0
dP=4
BC = rbc.Fixed(U, dP)

n_Ns = 5
Ns = np.zeros(n_Ns)

errs_p4 = np.zeros(n_Ns)
errs_p2 = np.zeros(n_Ns)
dN = 2
N = 20

# analytic solution dP for LambdaBump
fac_lam = np.sqrt(1-lam)

analytic_pert2_dP = -12*((np.pi*lam)**2)/(5*fac_lam**3)

analytic_pert4_dP = -8*(np.pi**4)* (428*(-1+fac_lam) - 214*(-2 + fac_lam)*lam - 53*(lam**2))/(175*fac_lam)


 
for k in range(n_Ns):
    Ns[k] = N 
    
    solver = control.Reynolds_Solver(Example, BC, args)
    pert = solver.fd_pert_solve(N,order=4, plot=False, get_all=True)

    
    if pert.order >1:

        errs_p2[k] = abs(pert.dP_pert2 -analytic_pert2_dP)/abs(analytic_pert2_dP)

    if pert.order > 2:
        # print(pert.dP_pert4)

        errs_p4[k] = abs(pert.dP_pert4 -analytic_pert4_dP)/abs(analytic_pert4_dP)
    
    N *= dN
    
linthresh = 1e-8
graphics.plot_log_multi([errs_p2, errs_p4], Ns, f'ELT convergence $(\lambda={lam}, \delta={delta})$', 
                        ['$\Delta P_2$', '$\Delta P_4$'], ['N', '$\Delta P_k$ Rel. Error '], linthresh, bigO_on=True, loc='left')


p4_rate = convergence.convg_rate(errs_p4)
p2_rate = convergence.convg_rate(errs_p2)
cnvg_rates = np.stack([p2_rate, p4_rate], axis=0)


    
print("cnvg rates")
print("p4: " + np.array2string(p4_rate))
print("p2: " + np.array2string(p2_rate))