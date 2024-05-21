#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Nov 12 12:53:47 2022

Graphics helpers
@author: sarahdennis
"""
import numpy as np
from matplotlib import pyplot as pp
from matplotlib import colors
from matplotlib import patches

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
    pp.plot(xs, fs, color='b', linewidth=.8)

    pp.title(title, fontweight="bold")
    
    pp.xlabel(axis[0])
    
    pp.ylabel(axis[1])
    # pp.ylim(0, 1.1*max(fs))
    
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


def plot_stream(vx, vy, xs, ys, title, ax_labels):
    
    pp.rcParams['figure.dpi'] = 500
    pp.figure()
    
    X, Y = np.meshgrid(xs, ys)
    
    m = len(ys)/len(xs)
    h_max = max(ys)
    stream_density_unbroken=[1,.6*m] 
    
    pp.streamplot(xs, ys, vx, vy, stream_density_unbroken, linewidth=0.5, color='k', broken_streamlines=False)
    
    #remove arrows
    ax = pp.gca()
    for art in ax.get_children():
        if not isinstance(art, patches.FancyArrowPatch):
            continue
        art.remove()        
    
    pp.title(title, fontweight="bold")
    pp.xlabel(ax_labels[0])
    pp.ylabel(ax_labels[1])

    ax.set_aspect('equal')
    ax.set_ylim(0,1.01*h_max)
    pp.show()
    
        
def plot_stream_heat(vx, vy, xs, ys, psi, title, ax_labels):
    
    pp.rcParams['figure.dpi'] = 500
    pp.figure()

    X, Y = np.meshgrid(xs, ys)
    
    stream_density=[1,2] #len(ys) = 2 len(xs)
    
    norm_symLog = colors.SymLogNorm(linthresh=1e-12, linscale=0.35, vmin=-1, vmax=1, clip=True)

    stream_plot=pp.streamplot(xs, ys, vx, vy, stream_density, broken_streamlines=False, linewidth=0.5, color=psi, cmap='Spectral_r', norm=norm_symLog)
    pp.colorbar(stream_plot.lines, label=ax_labels[0])
    ax = pp.gca()
    for art in ax.get_children():
        if not isinstance(art, patches.FancyArrowPatch):
            continue
        art.remove()        
    

    pp.title(title, fontweight="bold")
    pp.xlabel(ax_labels[1])
    pp.ylabel(ax_labels[2])

    ax.set_aspect('equal')
    ax.set_ylim(0,2) #min max ys
    pp.show()

def plot_stream_height(vx, vy, hs, xs, ys, title, ax_labels):
    
    pp.rcParams['figure.dpi'] = 500
    pp.figure()
    
    X, Y = np.meshgrid(xs, ys)
    
    m = len(ys)/len(xs)
    h_max = max(ys)
    stream_density_unbroken=[1,.6*m] 
    
    
    pp.streamplot(xs, ys, vx, vy, stream_density_unbroken, linewidth=0.5, color='k', broken_streamlines=False)
    
    ax = pp.gca()
    for art in ax.get_children():
        if not isinstance(art, patches.FancyArrowPatch):
            continue
        art.remove()        
    
    pp.plot(xs, hs, linewidth=0.8, color='r', label='$h(x)$')
    
    pp.title(title, fontweight="bold")
    pp.xlabel(ax_labels[0])
    pp.ylabel(ax_labels[1])

    ax.set_aspect('equal')
    ax.set_ylim(0,h_max)
    # pp.legend(loc='upper left')
    pp.show()
    
def plot_quiver_height(vx, vy, hs, xs, ys, title, ax_labels):
    
    pp.rcParams['figure.dpi'] = 500
    pp.figure()
    
    X, Y = np.meshgrid(xs, ys)
    
    m = len(ys)
    n = len(xs)
    h_max = max(ys)
    
    quiver_density = 50
    vx = mask(vx, ys, hs, m, n, quiver_density)
    vy = mask(vy, ys, hs, m, n, quiver_density)
    
    pp.quiver(xs, ys, vx, vy, color='k', scale=8, scale_units='x')
    pp.plot(xs, hs, linewidth=0.8, color='r', label='$h(x)$')
    
    pp.title(title, fontweight="bold")
    pp.xlabel(ax_labels[0])
    pp.ylabel(ax_labels[1])
    ax = pp.gca()

    # ax.set_aspect('equal')
    ax.set_ylim(0,1.25*h_max)
    pp.legend(loc='upper left')
    pp.show()

def mask(grid, ys, hs, m, n, density):
    mask = grid.copy()
    
    for j in range(m):
        for i in range(n):
            if ys[j] > hs[i]:
                mask[j,i] = None
            elif i % density != 0 or j%density != 0:
                    mask[j,i] = None
    return mask
        

def plot_contour(zs, xs, ys, title, labels):
    pp.rcParams['figure.dpi'] = 500
    pp.figure()

    X, Y = np.meshgrid(xs, ys)
    n_contours = max(zs.shape)
    
    contour_plot = pp.contour(X, Y, zs, n_contours, cmap='Spectral_r')
    
    pp.title(title, fontweight="bold")
    pp.xlabel(labels[1])
    pp.ylabel(labels[2])
    pp.colorbar(contour_plot, label=labels[0])
    
    ax = pp.gca()
    ax.set_aspect('equal', 'box')
    pp.show()

    
def plot_heatMap(zs, xs, ys, title, labels):
    pp.rcParams['figure.dpi'] = 500
    pp.figure()
    
    X, Y = np.meshgrid(xs, ys)
    color_plot = pp.pcolor(X, Y, zs, cmap='Spectral_r', norm=colors.SymLogNorm(linthresh=1e-8, linscale=0.25))
    pp.colorbar(color_plot, label=labels[0])

    pp.title(title, fontweight="bold")
    pp.xlabel(labels[1])
    pp.ylabel(labels[2])
    ax = pp.gca()

    ax.set_aspect('equal', 'box')
    pp.show()

def plot_contour_heat(zs, xs, ys, title, labels):
    pp.rcParams['figure.dpi'] = 500
    pp.figure()
    
    X, Y = np.meshgrid(xs, ys)

    norm_symLog = colors.SymLogNorm(linthresh=1e-12, linscale=0.35, vmin=-1, vmax=1, clip=True)

    color_plot = pp.pcolor(X, Y, zs, cmap='Spectral_r', norm=norm_symLog)
    
    pp.colorbar(color_plot, label=labels[0])
    

    n_contours = 10

    pp.rcParams["lines.linewidth"] = .15
    pp.contour(X, Y, zs, n_contours, colors='white')
    
    pp.title(title, fontweight="bold")
    pp.xlabel(labels[1])
    pp.ylabel(labels[2])
    
    ax = pp.gca()
    ax.set_aspect('equal', 'box')
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
    fig.legend(bbox_to_anchor=(1.1, 0.7))
    
    return fig

