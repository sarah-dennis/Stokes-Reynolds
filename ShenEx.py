# -*- coding: utf-8 -*-
"""
Created on Mon Oct 23 13:40:16 2023

 Comparison with Shen et al. 2018

@author: sarah
"""
import domain as dm
import examples as ex
import velocity as vel
import graphics as graph
#boundary velocity
U = 1 # inlet 

#boundary pressures
p0 = 0
pN = 0

#viscosity (dynamic)
visc = 0.188 #Pa*s 

#y=h(x) params
h_max = 250e-3 #mm
h_min = 133.976e-3 #mm

#x params
l_in = int(8.975e3) #mme3
l_tot = int(12.5e3) #mme3
x0 = 0
x_step = x0 + l_in
xf = x0 + l_tot

BC = 'fixed'
Nx = 500

domain = dm.Domain(x0, xf, visc, U, Nx, BC)

height, pressure = ex.step(domain, p0, pN, x_step, h_max, h_min)
# height, pressure = ex.variableStepLen_schurLUSolve(domain, p0, pN, x_step, h_max, h_min)

p_h_title = "Pressure (%s) and Height for %s"%(pressure.p_str, height.h_str)
p_h_labels = ["Pressure $p(x)$", "Height $h(x)$", "$x$"]
graph.plot_2D_twin(pressure.ps, height.hs, domain.xs, p_h_title, p_h_labels)



velocity = vel.Velocity(domain, height, pressure)
phv_title = "Pressure and Velocity for %s"%height.h_str
phv_fun_labels = ['velocity $(v_x, v_y)$', 'pressure $p(x)$', 'height $h(x)$']
phv_ax_labels =  ['$x$', '$y$']

graph.plot_phv(pressure.ps, height.hs, velocity.vx, velocity.vy, domain.xs, domain.ys, phv_title, phv_fun_labels,  phv_ax_labels)


print("max pressure: %.3f"%pressure.ps[domain.get_index(x_step)])
# expect max pressure 54 kpa
# reports 53967.503 (e-3 kpa)