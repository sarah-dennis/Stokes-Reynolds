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
    pp.title(title)
    fig.legend()
    return fig
    
def plot_2D(fs, xs, title, label):
    fig = pp.figure()
    pp.plot(xs, fs, label=label, color='g')
    
    
    #ax.set_xlim([0, 1])
    #ax.set_ylim([0, 1])

    pp.xlabel('x')
    #ax.set_ylabel()
    pp.title(title)
    fig.legend()
    return fig