# -*- coding: utf-8 -*-
"""
Created on Tue Oct 31 15:07:45 2023

BFS hydraulic resistance
Stokes vs Reynolds comparison

@author: sarah
"""
import graphics as graph
import domain as dm
from heights import StepHeight
from pressures import StepPressure
import numpy as np
# x boundary
x0 = 0
xf = 4
BC = "fixed" 

#pressure boundary
p0 = 2
pN = 1

#fluid properties
visc = 1 # kinematic viscosity
rho = 1 #density
lam = 48 #friction ??

# Velocity boundary
U = 0 #V_x(x,0) lower surface 

# Domain
Nx = 500
domain = dm.Domain(x0, xf, visc, U, Nx, BC)

#Height
l1 = 1
x_step = x0+l1
l2 = xf-l1

h1 = 1
H_ratio = 1.9423 #expansion ratio H/h 
h2 = H_ratio * h1

h_avg = (h1 + h2)/2
r = abs(h1 - h2)/2
height = StepHeight(domain, x_step, h1, h2)

#----------------------------------------------
# Reynolds analytic solution

p_reyn = StepPressure(domain, height, p0, pN)

px1_reyn = p_reyn.m_in
px2_reyn = p_reyn.m_out

p_max_reyn = px1_reyn*l1+p0
print("Reynolds \n px1: %.3f, px2: %.3f"%(px1_reyn, px2_reyn))
print(" p_max: %.3f\n"%p_max_reyn)

#----------------------------------------------

# Flux 
Q = (U * h1)/2 - (px1_reyn * (h1**3))/(12*visc) 
print("Flux: Q=%.3f"%Q)

# Bulk velocity
Ub1 = Q/h1
Ub2 = Q/h2 # = h1/h2 * Ub1

#reynolds number
Re = 2*h1*rho*Ub1/visc 
print("Reynolds num: %.3f\n"%Re)

#----------------------------------------------
# BISWAS (stokes) solution

dp1_biswas = -lam*rho*l1*(Ub1**2)/(2*h1*Re)   
dp2_biswas = -lam*rho*l2*(Ub2**2)/(2*h2*Re)

px1_biswas = dp1_biswas/l1
px2_biswas = dp2_biswas/l2

p_max_biswas = dp1_biswas+p0

# pressure drop adjusted for expansion ratio
px_corr = (1-h1**2/h2**2)*(rho*Ub1**2)/2

print("Biswas Stokes \n px1: %.3f, px2: %.3f "%(px1_biswas, px2_biswas))

print(" px_correction: %.3f\n"%px_corr)

print("hyd-res reyn: %.5f"%(-Q/(px1_reyn + px2_reyn)))
print("hyd-res biswas: %.5f"%(-Q/(px1_biswas + px2_biswas - px_corr)))

def make_disc_ps(x0, l1, l2, m1, m2, Nx):
    ps = np.zeros(Nx)
    dx = (l1 + l2)/Nx
    
    for i in range(Nx):
        xi = x0 + i*dx
        
        if xi < l1:
            ps[i] = p0 + i*m1*dx
        else:
            
            ps[i] = pN - (Nx-i)*m2*dx

    return ps

ps = make_disc_ps(x0, l1, l2, px1_biswas, px2_biswas, domain.Nx)

#------------------------------------------------------------------------------
# Pressure & Height plotting 
#------------------------------------------------------------------------------
p_h_title = "Reynolds Pressure and Height for BFS"
p_h_title2 = "Biswas Pressure and Height for BFS"
p_h_labels = ["Pressure $p(x)$", "Height $h(x)$", "$x$"]
graph.plot_2D_twin(p_reyn.ps, height.hs, domain.xs, p_h_title, p_h_labels)

graph.plot_2D_twin(ps, height.hs, domain.xs, p_h_title2, p_h_labels)





