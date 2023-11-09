# -*- coding: utf-8 -*-
"""
Created on Tue Oct 31 15:07:45 2023

BFS hydraulic resistance
Stokes vs Reynolds comparison

@author: sarah
"""
import _graphics as graph
import domain as dm
from heights import StepHeight
from pressures import StepPressure
import numpy as np
# x boundary
x0 = 0
xf = 4
BC = "fixed" 

#pressure boundary
p0 = 20
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
l1 = 0.5
x_step = x0+l1
l2 = xf-l1

h1 = 1
H_ratio = 1.9423 #expansion ratio 
h2 = H_ratio * h1

h_avg = (h1 + h2)/2
r = abs(h1 - h2)/2
height = StepHeight(domain, x_step, h_avg, r)

#reynolds analytic solution
p_reyn = StepPressure(domain, height, p0, pN)

px1_reyn = p_reyn.m_in
px2_reyn = p_reyn.m_out
p_max_reyn = px1_reyn*l1+p0
print("Reynolds \n px1: %.3f, px2: %.3f"%(px1_reyn, px2_reyn))
print(" p_max: %.3f\n"%p_max_reyn)

# Flux 
Q = (U * h1)/2 - (px1_reyn * (h1**3))/(12*visc) 
print("Flux: Q=%.3f"%Q)
# Bulk velocity
Ub1 = Q/h1
Ub2 = Q/h2 # = h1/h2 * Ub1

# BISWAS (stokes) solution

Re = 2*h1*rho*Ub1/visc #reynolds number

print("Reynolds num: %.3f\n"%Re)

px1_stokes = -lam*rho*l1*(Ub1**2)/(2*h1*Re*l1)   
p_max_stokes_in = px1_stokes*l1+p0


px2_stokes = -lam*rho*l2*(Ub2**2)/(2*h2*Re*l2)
p_max_stokes_out = -px2_stokes*(l2) + pN


#adjust for expansion ratio
px_corr = (1-h1**2/h2**2)*(rho*Ub1**2)/2

print("Biswas Stokes \n px1: %.3f, px2: %.3f "%(px1_stokes, px2_stokes))
print(" p_max_1, p_max_2: %.3f, %.3f\n"%(p_max_stokes_in, p_max_stokes_out))
print("px_correction: %.3f\n"%px_corr)

print("hyd-res reyn: %.5f"%(-Q/(px1_reyn + px2_reyn)))
print("hyd-res stokes: %.5f"%(-Q/(px1_stokes + px2_stokes - px_corr)))

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

ps = make_disc_ps(x0, l1, l2, px1_stokes, px2_stokes, domain.Nx) #p_reyn.ps

#------------------------------------------------------------------------------
# Pressure & Height plotting 
#------------------------------------------------------------------------------
p_h_title = "Reynolds Pressure and Height for BFS"
p_h_title2 = "Biswas Pressure and Height for BFS"
p_h_labels = ["Pressure $p(x)$", "Height $h(x)$", "$x$"]
graph.plot_2D_twin(p_reyn.ps, height.hs, domain.xs, p_h_title, p_h_labels)

graph.plot_2D_twin(ps, height.hs, domain.xs, p_h_title2, p_h_labels)





