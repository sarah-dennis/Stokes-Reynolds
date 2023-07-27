#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Nov 12 12:53:47 2022

Graphics helpers
@author: sarahdennis
"""
import numpy as np
from matplotlib import pyplot as pp
import csv
    
theta = 30
phi = 30


def plot_3D(f_2D, xs, zs, title):
    
    X, Z = np.meshgrid(xs, zs)
    
    pp.figure()
    ax = pp.axes(projection='3d')
    ax.plot_surface(X.T, Z.T, f_2D, rstride=1, cstride=1,cmap='viridis')
    pp.title(title)
    pp.xlabel('x')
    pp.ylabel('z')
    ax.view_init(theta, phi)

def plot_2D_multi(fs, xs, title, labels, axis):
    fig = pp.figure()
    ax = fig.add_subplot()
    colors = ['r', 'b', 'g', 'o', 'p']

    for i in range(len(fs)):
        
        ax.plot(xs, fs[i], label=labels[i], color=colors[i])
    
    #ax.set_xlim([0, 1])
    #ax.set_ylim([0, 1])

    ax.set_xlabel(axis[0])
    ax.set_ylabel(axis[1])
    pp.title(title,  fontweight ="bold")
    fig.legend()
    return fig
    
def plot_2D(fs, xs, title, y_label, x_label):
    fig = pp.figure()
    pp.plot(xs, fs, label=y_label, color='b')

    pp.title(title)
    
    pp.xlabel(x_label)
    
    pp.ylabel(y_label)
    pp.ylim(0, 1.1*max(fs))
    #fig.legend()
    
    return fig

def plot_2D_twin(fs, gs, xs, title, labels):
    pp.rcParams['figure.dpi'] = 300
    # Creating plot with dataset_1
    fig, ax1 = pp.subplots()
    
    ax1.set_xlabel(labels[2])

    color = 'tab:red'
    ax1.plot(xs, fs, color = color)
    ax1.set_ylabel(labels[0], color = color)
    ax1.tick_params(axis ='y', labelcolor = color)

    # Adding Twin Axes to plot dataset_2
    ax2 = ax1.twinx()
     
    color = 'tab:blue'
    ax2.set_ylabel(labels[1], color = color)
    ax2.plot(xs, gs, color = color)
    ax2.tick_params(axis ='y', labelcolor = color)
    ax2.set_ylim(0, 1.2*max(gs))
     
    # Adding title
    pp.title(title, fontweight ="bold")
     
    # Show plot
    pp.show()
    

def plot_log(fs, xs, title, y_label, x_label):
    fig = pp.figure()
    pp.loglog(xs, fs, color='b')   

    pp.title(title)
    
    pp.xlabel(x_label)
    
    pp.ylabel(y_label)
    
    return fig

    
def plot_log_multi(fs, xs, title, f_labels, axis):
    pp.rcParams['figure.dpi'] = 300
    # pp.rcParams["legend.loc"] = 'center right'
    fig = pp.figure()
    
    
    ax = fig.add_subplot()
    colors = ['r', 'b', 'g', 'purple']
    for i in range(len(fs)):
        
        ax.loglog(xs, fs[i], label=f_labels[i], color=colors[i])
    
    #ax.set_xlim([0, 1])
    #ax.set_ylim([0, 1])

    ax.set_xlabel(axis[0])
    ax.set_ylabel(axis[1])
    pp.title(title,  fontweight ="bold")
    fig.legend( bbox_to_anchor=(1.22, 0.5))
    return fig

