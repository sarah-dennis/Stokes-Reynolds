#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Nov 12 12:53:47 2022

Graphics helpers
@author: sarahdennis
"""
import numpy as np
from matplotlib import pyplot as pp

    
fig_theta = 10
fig_phi = 30

fig_Nx=100
fig_Nz=100


def plot_3D(f_2D, xs, zs, title):
    
    X, Y = np.meshgrid(xs, zs)
    
    pp.figure()
    ax = pp.axes(projection='3d')
    ax.plot_surface(X, Y, f_2D, rstride=1, cstride=1,cmap='viridis')
    pp.title(title)
    pp.xlabel('x')
    pp.ylabel('z')
    ax.view_init(fig_theta, fig_phi)
    
    
def plot_2d(f, xs):
    pp.figure()
    pp.plot(xs, f, label="$y=f(x)$")
    pp.title('title')
    pp.xlabel('x')
    pp.ylabel('y')
    pp.legend()
