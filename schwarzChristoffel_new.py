#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 30 11:05:44 2023

@author: sarahdennis
"""
import numpy as np
import cmath
from matplotlib import pyplot as pp

 
#------------------------------------------------------------------------------
# Pentagon BFS
#------------------------------------------------------------------------------
# -  M ----------- L -
# h1 |             | |
# -  N ---- P      | h2
#           |      | |
#           R ---- S -     
#    |--l1--|--l2--|  

# flow:  <--U-- (U>0)
U= 1/2

l1 = 1
l2 = 4

h1 = 1
h2 = 2

#-----------------------------------------------------------   
# # Schwarz-Christoffel Mapping z = f(w) where w = psi + i phi
#------------------------------------------------------------------------------       

a = (h2/h1)**2
b1 = h2/np.pi 
b2 = h1/np.pi
K = np.pi/(U*h2)
S = (U*h2)/np.pi
def f(w):
    eKw = np.exp(K*w.real)
    sinKw = cmath.sin(K*w.imag)
    cosKw = cmath.cos(K*w.imag)

    u1 = ( 2 * eKw   * cosKw - a-1)/(a-1)
    v1 = ( 2 * eKw   * sinKw      )/(a-1)
    u2 = (-2 * a/eKw * cosKw + a+1)/(a-1)
    v2 = ( 2 * a/eKw * sinKw      )/(a-1)
    
    alph1 = cmath.sqrt((1+u1)**2 + v1**2)
    beta1 = cmath.sqrt((1-u1)**2 + v1**2)
    alph2 = cmath.sqrt((1+u2)**2 + v2**2)
    beta2 = cmath.sqrt((1-u2)**2 + v2**2)

    acosh1_real = cmath.acosh((alph1 + beta1)/2)
    acosh1_imag = cmath.acosh((alph1 - beta1)/2) 
    acosh2_real = cmath.acosh((alph2 + beta2)/2)
    acosh2_imag = cmath.acosh((alph2 - beta2)/2)

    f_real = b1 * acosh1_real - b2 * acosh2_real
    f_imag = b1 * acosh1_imag - b2 * acosh2_imag - complex(0,(h2-h1))
    return f_real + f_imag

#------------------------------------------------------------------------------\
#sample w from psi-phi grid    
N_psi = 200
N_phi = 200


# psi_min = 0.441271 
# psi_max = 0.4412715
# phi_min = 0 
# phi_max = 0.0000005

psi_min = -1
psi_max = 2
phi_min = 0
phi_max = 1

psis = np.linspace(psi_min, psi_max, N_psi)
phis = np.linspace(phi_min, phi_max, N_phi)

#------------------------------------------------------------------------------
f_zs = np.zeros((N_psi, N_phi), complex)
for i in range(N_psi):
    for j in range(N_phi):

        w = complex(psis[i], phis[j])
        f_zs[i, j] = f(w)
        
#f_zs[i,:] = {z : psi(z) = psi[i]}
#f_zs[:,j] = {z : phi(z) = phi[j]}
# ------------------------------------------------------------------------------

print("Complete: starting plot")

#------------------------------------------------------------------------------
# Plotting
#------------------------------------------------------------------------------
title = "Potential $\psi$ & Stream $\phi$ "
subtitle = "$\psi(z) \in [%.3f, %.3f]$ \n $\phi(z) \in [%.1f, %.3f]$"%(psi_min, psi_max, phi_min, phi_max)
ax =["x","y"]
pp.figure()
pp.xlabel(ax[0])
pp.ylabel(ax[1])

pp.title(title + '\n' + subtitle)

for j in range(N_phi):
    phi_contr, = pp.plot(f_zs[:,j].real, f_zs[:,j].imag, color="indigo")
    
for i in range(N_psi):
    psi_contr, = pp.plot(f_zs[i].real, f_zs[i].imag, color="coral")

pp.legend( [phi_contr, psi_contr], ["$\phi$", "$\psi$"])

















