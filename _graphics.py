#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Nov 12 12:53:47 2022

Graphics helpers
@author: sarahdennis
"""
import numpy as np
from matplotlib import pyplot as pp

#------------------------------------------------------------------------------
def plot_3D(f_2D, xs, zs, title):
            
    theta = 30
    phi = 30
    X, Z = np.meshgrid(xs, zs)
    
    pp.figure()
    ax = pp.axes(projection='3d')
    ax.plot_surface(X.T, Z.T, f_2D, rstride=1, cstride=1,cmap='viridis')
    pp.title(title)
    pp.xlabel('x')
    pp.ylabel('z')
    ax.view_init(theta, phi)
    
#------------------------------------------------------------------------------
def plot_2D(fs, xs, title, axis):
    fig = pp.figure()
    pp.plot(xs, fs, color='b')

    pp.title(title)
    
    pp.xlabel(axis[0])
    
    pp.ylabel(axis[1])
    pp.ylim(0, 1.1*max(fs))
    
    return fig

def scatter_2D(fs, xs, title, axis):
    fig = pp.figure()
    pp.scatter(xs, fs, color='b')

    pp.title(title)
    
    pp.xlabel(axis[0])
    
    pp.ylabel(axis[1])
    pp.ylim(0, 1.1*max(fs))
    
    return fig


def plot_2D_multi(fs, xs, title, fun_labels, ax_labels):
    fig = pp.figure()
    ax = fig.add_subplot()
    colors = ['r', 'b', 'g', 'orange', 'purple']

    for i in range(len(fs)):
        
        ax.plot(xs, fs[i], label=fun_labels[i], color=colors[i])
    
    #ax.set_xlim([0, 1])
    #ax.set_ylim([0, 1])

    ax.set_xlabel(ax_labels[0])
    ax.set_ylabel(ax_labels[1])
    pp.title(title,  fontweight ="bold")
    fig.legend()
    return fig

#------------------------------------------------------------------------------
def plot_2D_twin(fs, gs, xs, title, ax_labels):
    pp.rcParams['figure.dpi'] = 300
    # Creating plot with dataset_1
    fig, ax1 = pp.subplots()
    
    ax1.set_xlabel(ax_labels[2])
    color = 'tab:red'
    ax1.plot(xs, fs, color = color)
    ax1.set_ylabel(ax_labels[0], color = color)
    ax1.tick_params(axis ='y', labelcolor = color)

    # Adding Twin Axes to plot dataset_2
    ax2 = ax1.twinx()
    color = 'tab:blue'
    ax2.set_ylabel(ax_labels[1], color = color)
    ax2.plot(xs, gs, color = color)
    ax2.tick_params(axis ='y', labelcolor = color)
    ax2.set_ylim(0, 1.2*max(gs))
     
    pp.title(title, fontweight ="bold")
     
    pp.show()
    
#------------------------------------------------------------------------------
def plot_phv(ps, hs, vx, vy, xs, ys, title, fun_labels, ax_labels):   
    pp.rcParams['figure.dpi'] = 500
    
    fig = pp.figure()
    X, Y = np.meshgrid(xs, ys)
    
    # mask y > h(x)
    mask = np.zeros((len(xs), len(ys)), dtype=bool)
    for i in range(len(ys)):
        for j in range(len(xs)):
            mask[i,j] = ys[i] > hs[j]
    
    #Spectral pressure plot 
    P = [ps for y in ys]
    P = np.ma.array(P, mask=mask)
        
    press_plot = pp.scatter(X, Y, c=P, cmap=pp.cm.get_cmap('Spectral'))
    fig.colorbar(press_plot, label="pressure")

    #Velocity vector plot
    skip = 20
    
    thin_xs = xs[::skip]
    thin_ys = ys[::skip]
    
    thin_X, thin_Y = np.meshgrid(thin_xs, thin_ys)

    thin_vx = np.zeros((len(thin_xs), len(thin_ys)))
    thin_vy = np.zeros((len(thin_xs), len(thin_ys)))
    thin_mask = np.zeros((len(thin_xs), len(thin_ys)))
    j = 0
    for i in range(0, len(vx), skip):
        thin_vx[j] = vx[i][::skip]
        thin_vy[j] = vy[i][::skip]
        thin_mask[j] = mask[i][::skip]
        j+=1
    
    thin_vx = np.ma.array(thin_vx, mask=thin_mask)
    thin_vy = np.ma.array(thin_vy, mask=thin_mask)
    
    pp.quiver(thin_X, thin_Y, thin_vx, thin_vy, width=0.003)
    
    # pp.plot(xs, hs, label='height', color='white')
    
    # pp.legend()
    pp.title(title, fontweight="bold")
    pp.xlabel(ax_labels[0])
    pp.ylabel(ax_labels[1])
    

    pp.show()
    
#------------------------------------------------------------------------------   

def plot_log(fs, xs, title, ax_labels):
    fig = pp.figure()
    pp.loglog(xs, fs, color='b')   

    pp.title(title)
    
    pp.xlabel(ax_labels[0])
    pp.ylabel(ax_labels[1])
    
    return fig

    
def plot_log_multi(fs, xs, title, f_labels, ax_labels):
    pp.rcParams['figure.dpi'] = 300
    fig = pp.figure()
    
    
    ax = fig.add_subplot()
    colors = ['r', 'b', 'g', 'orange', 'purple']
    
    for i in range(len(fs)):
        ax.loglog(xs, fs[i], label=f_labels[i], color=colors[i])
        
    ax.set_xlabel(ax_labels[0])
    ax.set_ylabel(ax_labels[1])
    
    pp.title(title,  fontweight ="bold")
    fig.legend(bbox_to_anchor=(1.22, 0.5))
    
    return fig

