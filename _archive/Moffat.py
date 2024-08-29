# -*- coding: utf-8 -*-
"""
Created on Tue Nov 14 15:43:27 2023

@author: sarah
"""
import numpy as np
import domain 
import heights 
import velocity
import pressures
import _graphics 
from matplotlib import pyplot as pp

x0 = 0
xf = 5

p0 = 0
pN = 0

U = 5

visc = 1

Nx = 80
Ny = 100
BC = "fixed"

domain = domain.Domain(x0, xf, visc, U, Nx, BC)
 
alpha = np.pi/4


#height
x_step = 1
h1 = 0.5
h2 = 1
height = heights.StepHeight(domain, x_step, h1, h2)
domain.set_ys(height, Ny) 

# reynolds pressure...
pressure = pressures.StepPressure(domain, height, p0, pN)

p_h_title = "Pressure (%s) and Height for %s"%(pressure.p_str, height.h_str)
p_h_labels = ["Pressure $p(x)$", "Height $h(x)$", "$x$"]
_graphics.plot_2D_twin(pressure.ps, height.hs, domain.xs, p_h_title, p_h_labels)

# XY_polar = [y, x, (r, theta)]
XY_polar_grid = np.zeros((Ny, Nx, 2))

for j in range(Ny):
    y = domain.ys[j] #-h2
    
    for i in range(Nx):    
        x = domain.xs[i]#-x_step
        
        #r = sqrt(x^2 + y^2)
        XY_polar_grid[j,i][0] = (x**2 + y**2)**(1/2)
        
        #theta = arctan(y/x)
        if x == 0:
            XY_polar_grid[j,i][1] = np.pi/2
        else:
            XY_polar_grid[j,i][1] = np.arctan(y/x) 


velx = np.zeros((Ny, Nx))
vely = np.zeros((Ny, Nx))


c = 1/(1-alpha**2)

for j in range(Ny):
    for i in range(Nx):
        r = XY_polar_grid[j,i][0]
        t = XY_polar_grid[j,i][1]
        
        sint = np.sin(t)
        cost = np.cos(t)
      
        f = (1-c)*sint + c*t*cost + alpha*c*t*sint
        df = (1-c)*cost+ c*(cost - t*sint) + alpha*c*(sint + t*cost)
        
        # if domain.xs[i] < x_step and domain.ys[j] < h1:
        velx[j, i] = U*df
        vely[j, i] = -U*f
        
        # elif domain.xs[i] > x_step and domain.ys[j]<h2 :
            # velx[j, i] = U*df
            # vely[j, i] = -U*f
            

vel_title = "Moffat Velocity in a corner alpha = %.2f"%alpha
vel_ax_labels = ['x', 'y']
_graphics.plot_vel(velx, vely, domain.xs, domain.ys, vel_title, vel_ax_labels)

reynVelocity = velocity.Velocity(domain, height, pressure)

phv_title = "Pressure and Velocity for %s"%height.h_str
phv_fun_labels = ['velocity $(v_x, v_y)$', 'pressure $p(x)$', 'height $h(x)$']
phv_ax_labels =  ['$x$', '$y$']

_graphics.plot_phv(pressure.ps, height.hs, reynVelocity.vx, reynVelocity.vy, domain.xs, domain.ys, phv_title, phv_fun_labels,  phv_ax_labels)
