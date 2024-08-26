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
    ax.plot_surface(X.T, Z.T, f_2D, rstride=1, cstride=1, cmap='viridis')
    pp.title(title)
    pp.xlabel('x')
    pp.ylabel('z')
    ax.view_init(theta, phi)
    
#------------------------------------------------------------------------------
def plot_2D(fs, xs, title, axis_labels, color='b'):
    fig = pp.figure()
    pp.plot(xs, fs, color='b', linewidth=.8)

    pp.title(title, fontweight="bold")
    
    pp.xlabel(axis_labels[0])
    
    pp.ylabel(axis_labels[1])
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
    colors = ['red', 'blue', 'orange', 'green', 'magenta']

    for i in range(len(fs)):
        
        ax.plot(xs, fs[i], label=fun_labels[i], color=colors[i])
    
    #ax.set_xlim([0, 1])
    #ax.set_ylim([0, 1])

    ax.set_xlabel(ax_labels[0])
    ax.set_ylabel(ax_labels[1])
    pp.title(title,  fontweight ="bold")
    fig.legend()
    return fig


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
    
        
def plot_stream_heat(vx, vy, xs, ys, psi, title, ax_labels, log_cmap, linthresh=1e-18, contour_density=1e-2):
    
    pp.rcParams['figure.dpi'] = 500
    pp.figure()

    X, Y = np.meshgrid(xs, ys)
    
    stream_density=[8*len(ys)/len(xs),1]
    
    if log_cmap:
    
        norm_symLog = colors.SymLogNorm(linthresh, linscale=0.25, vmin=-1, vmax=1, clip=True)
        stream_plot=pp.streamplot(xs, ys, vx, vy, stream_density, broken_streamlines=False, linewidth=0.5, color=psi, cmap='viridis', norm=norm_symLog)
    else:
        stream_plot=pp.streamplot(xs, ys, vx, vy, stream_density, broken_streamlines=False, linewidth=0.5, color=psi, cmap='viridis')

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
    ax.set_ylim(0,max(ys)) #min max ys
    pp.show()

def plot_stream_height(vx, vy, hs, xs, ys, title, ax_labels):
    
    pp.rcParams['figure.dpi'] = 800
    pp.figure()
    
    X, Y = np.meshgrid(xs, ys)
    
    # m = len(ys)/len(xs)
    h_max = max(ys)
    # stream_density_unbroken=[1,.6*m] 
    
    
    pp.streamplot(xs, ys, vx, vy, linewidth=0.5, color='k', broken_streamlines=False)
    
    ax = pp.gca()
    # for art in ax.get_children():
    #     if not isinstance(art, patches.FancyArrowPatch):
    #         continue
    #     art.remove()        
    
    pp.plot(xs, hs, linewidth=0.8, color='k', label='$h(x)$')
    
    pp.title(title, fontweight="bold")
    pp.xlabel(ax_labels[0])
    pp.ylabel(ax_labels[1])

    ax.set_aspect(xs[-1]/ys[-1])
    ax.set_ylim(0,1.1*h_max)
    # pp.legend(loc='upper left')
    pp.show()

#------------------------------------------------------------------------------
    
def plot_quiver_height(vx, vy, hs, xs, ys, title, ax_labels):

    pp.rcParams['figure.dpi'] = 500
    pp.figure()
    
    X, Y = np.meshgrid(xs, ys)
    
    m = len(ys)
    n = len(xs)
    h_max = max(ys)
    ly = max(ys)
    lx = max(xs)
    v_scale = 20*np.max(vx)/(lx)
    
    vx = quiver_mask(vx, m, n, ly, lx)
    vy = quiver_mask(vy, m, n, ly, lx)
    
    pp.quiver(xs, ys, vx, vy, color='k', width=.001, scale=v_scale, scale_units='x')
    pp.plot(xs, hs, linewidth=0.8, color='r', label='$h(x)$')
    
    pp.title(title, fontweight="bold")
    pp.xlabel(ax_labels[0])
    pp.ylabel(ax_labels[1])
    ax = pp.gca()

    ax.set_aspect('equal')
    ax.set_ylim(0,h_max)
    pp.legend(loc='upper left')
    pp.show()
    
def plot_quiver(vx, vy, xs, ys, title, ax_labels):

    pp.rcParams['figure.dpi'] = 500
    pp.figure()
    
    X, Y = np.meshgrid(xs, ys)
    
    m = len(ys)
    n = len(xs)
    ly = max(ys)
    lx = max(xs)
    v_scale = 20*np.max(vx)/(lx)

    vx = quiver_mask(vx, m, n, ly, lx)
    vy = quiver_mask(vy, m, n, ly, lx)

    pp.quiver(xs, ys, vx, vy, scale=v_scale, scale_units='x', width=.001, color='k')
    
    pp.title(title, fontweight="bold")
    pp.xlabel(ax_labels[0])
    pp.ylabel(ax_labels[1])
    ax = pp.gca()

    ax.set_aspect('equal')
    ax.set_ylim(0,ly)
    pp.show()

def quiver_mask(grid, m, n, ly, lx):
    density = (ly*4,lx)
    mask = grid.copy()
    j_mod, i_mod = density
    for j in range(m):
        for i in range(n):
            if i % i_mod != 0 or j%j_mod != 0:
                    mask[j,i] = None
    return mask
 
#------------------------------------------------------------------------------       

def plot_contour(zs, xs, ys, title, labels, log_cmap=False, linthresh=1e-18):
    pp.rcParams["lines.linewidth"] = .4
    pp.rcParams['figure.dpi'] = 1000

    pp.figure()
    
    X, Y = np.meshgrid(xs, ys)
    n_contours = max(zs.shape)//2

    if log_cmap:
        norm_symLog = colors.SymLogNorm(linthresh, linscale=0.25, vmin=-1, vmax=1, clip=True)
        contour_plot = pp.contour(X, Y, zs,  n_contours, cmap='plasma', norm=norm_symLog)
    else:
        contour_plot = pp.contour(X, Y, zs,  n_contours, cmap='plasma')
        
    # contour_plot = pp.contour(X, Y, zs, n_contours, cmap='plasma')
        
    pp.title(title, fontweight="bold")
    pp.xlabel(labels[1])
    pp.ylabel(labels[2])
    pp.colorbar(contour_plot, label=labels[0])
    
    ax = pp.gca()
    # ax.set_aspect('equal')
    pp.show()



def plot_contour_mesh(zs, xs, ys, title, labels, log_cmap=True, linthresh=1e-18, n_contours=20):
    pp.rcParams['figure.dpi'] = 1000
    pp.figure()
    # zs = np.ma.masked_where(zs == 0, zs)
    X, Y = np.meshgrid(xs, ys)
#'Spectral_r'
    if log_cmap:
        norm_symLog = colors.SymLogNorm(linthresh, linscale=0.25, vmin=-1, vmax=1, clip=True)
        color_plot = pp.pcolor(X, Y, zs, cmap='viridis', norm=norm_symLog)
    else:
        color_plot = pp.pcolor(X, Y, zs, cmap='viridis')
    
    pp.colorbar(color_plot, label=labels[0])
    

    pp.rcParams["lines.linewidth"] = .2
    pp.contour(X, Y, zs, n_contours, colors='white')
    
    pp.title(title, fontweight="bold")
    pp.xlabel(labels[1])
    pp.ylabel(labels[2])
    
    ax = pp.gca()
    ax.set_aspect('equal')
    pp.show()    
#------------------------------------------------------------------------------   

def plot_log(fs, xs, title, ax_labels):
    fig = pp.figure()
    pp.loglog(xs, fs, color='b')   

    pp.title(title)
    
    pp.xlabel(ax_labels[0])
    pp.ylabel(ax_labels[1])
    
    return fig

    
def plot_log_multi(fs, xs, title, f_labels, ax_labels, linthresh=1e-6, O1=1e-2, O1half=1e-3/2, O2=1e-3):
    pp.rcParams['figure.dpi'] = 300
    fig = pp.figure()
    
    
    ax = fig.add_subplot()
    colors = ['red', 'blue', 'orange', 'green', 'magenta']

    pp.rcParams["lines.linewidth"] = .8
    for i in range(len(fs)):
        ax.plot(xs, fs[i], label=f_labels[i], color=colors[i], marker='x', markevery=1)
    
    
    # reference lines
    ax.plot(xs, [O1*x**-1 for x in xs], label="$\mathcal{O}(%s^{-1})$"%ax_labels[0], color='darkgrey')    
    # ax.plot(xs,[O1half*x**-np.sqrt(2) for x in xs], label="$\mathcal{O}(%s^{-\sqrt{2}})$"%ax_labels[0], color='grey')
    ax.plot(xs, [O2*x**-2 for x in xs], label="$\mathcal{O}(%s^{-2})$"%ax_labels[0], color='lightblue')
    
    ax.set_xscale('log')
    ax.set_yscale('symlog', linthresh=linthresh)
    ax.set_ylim(0,1)
    ax.set_xlabel(ax_labels[0])
    ax.set_ylabel(ax_labels[1])
    
    pp.title(title,  fontweight ="bold")
    fig.legend(bbox_to_anchor=(0.5, 0.4))
    
    return fig

