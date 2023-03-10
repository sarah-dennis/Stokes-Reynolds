#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Nov 12 12:53:47 2022

Graphics helpers
@author: sarahdennis
"""
import numpy as np
from matplotlib import pyplot as pp
    
theta = 30
phi = 30

Nx=100
Nz=100

def plot_3D(f_2D, xs, zs, title):
    
    X, Z = np.meshgrid(xs, zs)
    
    pp.figure()
    ax = pp.axes(projection='3d')
    ax.plot_surface(X.T, Z.T, f_2D, rstride=1, cstride=1,cmap='viridis')
    pp.title(title)
    pp.xlabel('x')
    pp.ylabel('z')
    ax.view_init(theta, phi)

def plot_2D_multi(fs, xs, title, labels):
    fig = pp.figure()
    ax = fig.add_subplot()
    colors = ['r', 'b', 'g', 'o', 'p']
    for i in range(len(fs)):
        
        ax.plot(xs, fs[i], label=labels[i], color=colors[i])
    
    #ax.set_xlim([0, 1])
    #ax.set_ylim([0, 1])

    ax.set_xlabel('x')
    #ax.set_ylabel()
    pp.title(title,  fontweight ="bold")
    fig.legend()
    return fig
    
def plot_2D(fs, xs, title, y_label, x_label):
    fig = pp.figure()
    pp.plot(xs, fs, label=y_label, color='g')

    pp.title(title)
    
    pp.xlabel(x_label)
    
    pp.ylabel(y_label)
    #fig.legend()
    
    return fig

def plot_p_h(ps, ps_num, hs, xs, ex_title):

    # Creating plot with dataset_1
    fig, ax1 = pp.subplots()
    
    ax1.set_xlabel('$x$')
    
    color = 'tab:red'
    ax1.plot(xs, ps, color = color, label="exact")
    
    color ='tab:purple'
    ax1.plot(xs, ps_num, color = color, label="numerical")
    
    ax1.set_ylabel('Pressure $p(x)$', color = 'black')
    ax1.tick_params(axis ='y', labelcolor = 'black')
    #ax1.set_ylim(0, 0.5)
    
    fig.legend()
     
    # Adding Twin Axes to plot dataset_2
    ax2 = ax1.twinx()
     
    color = 'tab:blue'
    ax2.set_ylabel('Height $h(x)$', color = color)
    ax2.plot(xs,hs, color = color)
    ax2.tick_params(axis ='y', labelcolor = color)
    ax2.set_ylim(0, 3.5)
     
    # Adding title
    pp.title('Pressure for %s'%ex_title, fontweight ="bold")
     
    # Show plot
    pp.show()
    
    
    
    
    
    
    
    
    