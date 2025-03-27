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

# colour_map_stream = 'viridis' 
# colour_map_stream = 'Spectral_r' 
# colour_map_mesh = 'plasma'
# colour_map_mesh='PiYG'
# colour_map_mesh='RdYlBu_r'
# colour_map_mesh = 'Spectral_r'
# colour_map_stream = 'plasma'

# stream_cmap = pp.cm.plasma(np.arange(pp.cm.plasma.N))
stream_cmap = pp.cm.Spectral_r(np.arange(pp.cm.Spectral_r.N))
stream_cmap[:,0:3] *= 1
colour_map_stream = colors.ListedColormap(stream_cmap)


mesh_cmap = pp.cm.RdYlBu_r(np.arange(pp.cm.RdYlBu_r.N))
mesh_cmap[:,0:3] *= 0.95
colour_map_mesh = colors.ListedColormap(mesh_cmap)

# colour_bar_scale=0.015 # very long figures, H=1.25, L=4
# colour_bar_scale=0.024 # long figures like H=2, L=4
colour_bar_scale=0.04 # almost square figures like H=2.75, L=4 or zooms

dpi=400

contour_width = 0.25
stream_width = 1
line_width = 1.5

SMALL_SIZE = 10
MEDIUM_SIZE = 12
BIGGER_SIZE = 16

pp.rc('font', size=SMALL_SIZE)          # controls default text sizes
pp.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
pp.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
pp.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
pp.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
pp.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
pp.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

linthresh = 1e-4
#------------------------------------------------------------------------------
def plot_2D(fs, xs, title, axis_labels, color='darkmagenta'):
    fig = pp.figure()
    pp.rcParams['figure.dpi'] = dpi
    # marker = None
    marker = 'o'
    pp.plot(xs, fs, color=color, linewidth=.8, marker=marker)

    pp.title(title, fontweight="bold")
    
    pp.xlabel(axis_labels[0])
    
    pp.ylabel(axis_labels[1])
    # pp.ylim(0, 1.1*max(fs))
    pp.minorticks_on()
    return fig

def plot_2D_multi(fs, xs, title, fun_labels, ax_labels, loc='upper', colors='pri'):
    fig = pp.figure()
    pp.rcParams['figure.dpi'] = dpi
    ax = fig.add_subplot()
    if colors== 'pri':
        cs = ['r','forestgreen','b', 'darkmagenta', 'darkorgange']
    else: 
        cs=['forestgreen', 'darkmagenta', 'darkorgange']
    markers = ['D', 'o', 's', '*', 'H', 'X']
    for i in range(len(fs)):
        
        ax.plot(xs, fs[i], label=fun_labels[i], color=cs[i], linewidth=0.8, marker=markers[i])
    
    #ax.set_xlim([0, 1])
    # ax.set_ylim([-1, 1])

    ax.set_xlabel(ax_labels[0])
    ax.set_ylabel(ax_labels[1])
    pp.title(title,  fontweight ="bold")
    pp.minorticks_on()
    if loc== 'upper':
        fig.legend(bbox_to_anchor=(0.9, 0.875))
    elif loc=='left':
        fig.legend(bbox_to_anchor=(0.4, 0.875))
    else:
        fig.legend(bbox_to_anchor=(0.9, 0.275))
    return fig

def plot_2D_multi_multi(fs, xs, title, fun_labels, ax_labels, loc, colors):
    fig = pp.figure()
    pp.rcParams['figure.dpi'] = dpi
    ax = fig.add_subplot()
    if colors=='sec':
        cs = ['forestgreen', 'darkmagenta', 'orange', 'firebrick', 'royalblue', 'crimson']
    else:
        cs = ['firebrick', 'mediumblue','crimson','royalblue','indianred',  'dodgerblue']
    markers = ['D', 'o', 's', '^', 'd', 'h']
    for i in range(len(fs)):
        
        ax.plot(xs[i], fs[i], label=fun_labels[i], color=cs[i], linewidth=0.8, marker=markers[i])
    
    #ax.set_xlim([0, 1])
    #ax.set_ylim([0, 1])

    ax.set_xlabel(ax_labels[0])
    ax.set_ylabel(ax_labels[1])
    pp.title(title,  fontweight ="bold")
    if loc== 'upper':
        
        fig.legend(bbox_to_anchor=(0.9, 0.875))
    else:
        fig.legend(bbox_to_anchor=(0.9, 0.275))
    pp.minorticks_on()
    return fig

#------------------------------------------------------------------------------
def plot_stream(vx, vy, xs, ys, title, ax_labels):
    
    pp.rcParams['figure.dpi'] = dpi
    pp.figure()
    
    X, Y = np.meshgrid(xs, ys)

    stream_density=[1,1]
    pp.streamplot(xs, ys, vx, vy, stream_density, linewidth=0.5, color='k', broken_streamlines=False)
    
    #remove arrows
    ax = pp.gca()
    # for art in ax.get_children():
    #     if not isinstance(art, patches.FancyArrowPatch):
    #         continue
    #     art.remove()        
    
    pp.title(title, fontweight="bold")
    pp.xlabel(ax_labels[0])
    pp.ylabel(ax_labels[1])

    ax.set_aspect('equal')
    ax.set_ylim(0,1.01*max(ys))
    pp.show()
    
        
def plot_stream_heat(vx, vy, xs, ys, color_map, title, ax_labels, vmin, vmax, log_cmap=False, linthresh=linthresh):
    
    pp.rcParams['figure.dpi'] = dpi
    
    pp.figure()
    
    X, Y = np.meshgrid(xs, ys)
    
    stream_density=[xs.shape[0]/ys.shape[0],1]
    
    if log_cmap:
        norm_symLog = colors.AsinhNorm(linthresh, vmin=vmin, vmax=vmax, clip=False)
        stream_plot=pp.streamplot(xs, ys, vx, vy, stream_density, broken_streamlines=False, linewidth=stream_width, color=color_map, cmap=colour_map_stream, norm=norm_symLog)
    else:
        no_norm = colors.CenteredNorm(vcenter=vmin + vmax/2, halfrange=vmax/2, clip=True)
        stream_plot=pp.streamplot(xs, ys, vx, vy, stream_density, broken_streamlines=False, linewidth=stream_width, color=color_map, cmap=colour_map_stream, norm=no_norm)

    cb=pp.colorbar(stream_plot.lines, label=ax_labels[0], fraction=colour_bar_scale, pad=0.025)
    ticks = np.linspace(vmin, vmax, num=5)
    cb.set_ticks(ticks)

    # remove contours arrows
    ax = pp.gca()
    for art in ax.get_children():
        if not isinstance(art, patches.FancyArrowPatch):
            continue
        art.remove()        


    pp.title(title, fontweight="bold")
    pp.xlabel(ax_labels[1])
    pp.ylabel(ax_labels[2])
    ax.set_aspect('equal')
    pp.minorticks_on()
    pp.show()
    
    
       
def plot_quiver(vx, vy, xs, ys, color_map, title, ax_labels, vmin, vmax, linthresh=linthresh):
    
    pp.rcParams['figure.dpi'] = dpi
    
    pp.figure()
    
    X, Y = np.meshgrid(xs, ys)
    
    N = 10
    vscale=10
    
    # if log_cmap:
        # norm_symLog = colors.AsinhNorm(linthresh, vmin=vmin, vmax=vmax, clip=False)
    pp.quiver(xs[:: N], ys[:: N], vx[::N, ::N], vy[::N, ::N], scale=vscale)#, color=color_map, cmap=colour_map_stream, norm=norm_symLog)
    # else:
        # no_norm = colors.CenteredNorm(vcenter=vmin + vmax/2, halfrange=vmax/2)
        # quiver_plot=pp.quiver(xs, ys, vx, vy, scale=vscale)#, color=color_map, cmap=colour_map_stream, norm=no_norm)

    # cb=pp.colorbar(quiver_plot, label=ax_labels[0], fraction=colour_bar_scale, pad=0.025)
    # ticks = np.linspace(vmin, vmax, num=5)
    # cb.set_ticks(ticks)

    # remove contours arrows
    ax = pp.gca()
    # for art in ax.get_children():
    #     if not isinstance(art, patches.FancyArrowPatch):
    #         continue
    #     art.remove()        


    pp.title(title, fontweight="bold")
    pp.xlabel(ax_labels[1])
    pp.ylabel(ax_labels[2])
    ax.set_aspect('equal')
    pp.minorticks_on()
    pp.show()

#------------------------------------------------------------------------------       

def plot_contour(zs, xs, ys, title, labels, log_cmap=False, linthresh=linthresh):
    pp.rcParams["lines.linewidth"] = .5
    pp.rcParams['figure.dpi'] = dpi

    pp.figure()
    
    X, Y = np.meshgrid(xs, ys)
    n_contours = max(zs.shape)

    if log_cmap:
        norm_symLog = colors.AsinhNorm(linthresh)#, vmin=-1, vmax=1, clip=True)
        contour_plot = pp.contour(X, Y, zs,  n_contours, cmap='plasma', norm=norm_symLog)
    else:
        contour_plot = pp.contour(X, Y, zs,  n_contours, cmap='plasma')
        
        
    pp.title(title, fontweight="bold")
    pp.xlabel(labels[1])
    pp.ylabel(labels[2])
    pp.colorbar(contour_plot, label=labels[0])
    
    ax = pp.gca()
    ax.set_aspect('equal')
    pp.show()
    
def plot_contour_multi(funs, xs, ys, title, fun_labels, labels, y_lim=None):
    pp.rcParams["lines.linewidth"] = .5
    pp.rcParams['figure.dpi'] = dpi

    pp.figure()
    
    X, Y = np.meshgrid(xs, ys)
    n_contours = len(xs)//2
    colors = [ 'forestgreen', 'darkorchid', 'darkorange', 'royalblue', 'firebrick']  
    plots = []
    
    for i in range(len(funs)):          
        pp.contour(X, Y, funs[i], n_contours, colors=colors[i])
        plots.append(pp.plot(0,0,color=colors[i], label=fun_labels[i]))
        
        
    pp.title(title, fontweight="bold")
    pp.xlabel(labels[0])
    pp.ylabel(labels[1])    
    ax = pp.gca()
    pp.legend()
    # ax.set_aspect('equal')
    if y_lim is not None:
        ax.set_ylim(0, y_lim)
    pp.show()    



def plot_contour_mesh(zs, xs, ys, title, labels, vmin, vmax, log_cmap=False, linthresh=linthresh, n_contours=30, y_lim=None):
    pp.rcParams['figure.dpi'] = dpi
    pp.figure()
    X, Y = np.meshgrid(xs, ys)

    if log_cmap:
        norm_symLog = colors.AsinhNorm(linthresh, vmin=vmin, vmax=vmax, clip=True)
        color_plot = pp.pcolor(X, Y, zs, cmap=colour_map_mesh, norm=norm_symLog)
    else:
        color_plot = pp.pcolor(X, Y, zs, cmap=colour_map_mesh, vmin=vmin, vmax=vmax)
    
    cb = pp.colorbar(color_plot, label=labels[0], fraction=colour_bar_scale, pad=0.025)
    
    
    ticks = np.linspace(vmin, vmax, num=5)
    cb.set_ticks(ticks)

    pp.rcParams["lines.linewidth"] = contour_width
    pp.contour(X, Y, zs, n_contours, colors='black')
    
    pp.title(title, fontweight="bold")
    pp.xlabel(labels[1])
    pp.ylabel(labels[2])

    pp.minorticks_on()
    ax = pp.gca()
    ax.set_aspect('equal')    
    ax.set_ylim(min(ys),y_lim)
    # ax.set_facecolor('black')
    pp.show()    
#------------------------------------------------------------------------------   

def plot_log(fs, xs, title, ax_labels):
    fig = pp.figure()
    pp.rcParams['figure.dpi'] = dpi
    pp.loglog(xs, fs, color='b')   

    pp.title(title)
    
    pp.xlabel(ax_labels[0])
    pp.ylabel(ax_labels[1])
    
    return fig


def plot_log_multi(fs, xs, title, f_labels, ax_labels, linthresh=linthresh, O1=1,O2=1e1):
    pp.rcParams['figure.dpi'] = dpi
    fig = pp.figure()
    
    
    ax = fig.add_subplot()
    colors = [ 'forestgreen', 'darkorchid', 'darkorange', 'royalblue', 'firebrick']    
    markers = ['D', 'o', 's', '*', 'H', 'X']

    pp.rcParams["lines.linewidth"] =line_width
    for i in range(len(fs)):
        ax.plot(xs, fs[i], label=f_labels[i], color=colors[i], marker=markers[i], markevery=1)
    
    
    # reference lines
    ax.plot(xs, [O1*x**-1 for x in xs], label="$\mathcal{O}(%s^{-1})$"%ax_labels[0], color='darkgrey')    
    ax.plot(xs, [O2*x**-2 for x in xs], label="$\mathcal{O}(%s^{-2})$"%ax_labels[0], color='k')
    
    ax.set_xscale('log')
    ax.set_yscale('symlog', linthresh=linthresh)

    ax.set_ylim(0, np.max([O1,2*np.max(fs)]))

    ax.set_xlabel(ax_labels[0])
    ax.set_ylabel(ax_labels[1])
    
    pp.title(title,  fontweight ="bold")
    # fig.legend(bbox_to_anchor=(0.3, 0.325))
    fig.legend(bbox_to_anchor=(0.3, 0.425))    
    return fig
#------------------------------------------------------------------------------------
def grid_zoom_2D(grid, ex, x_start, x_stop, y_start, y_stop):
    i_0 = int((x_start - ex.x0)/ex.dx)
    i_f = int((x_stop - ex.x0)/ex.dx)
    j_0 = int((y_start - ex.y0)/ex.dy)
    j_f = int((y_stop - ex.y0)/ex.dy)
    return grid[j_0:j_f,i_0:i_f]

def grid_zoom_1D(grid_x, grid_y, ex, x_start, x_stop, y_start, y_stop):
    i_0 = int((x_start - ex.x0)/ex.dx)
    i_f = int((x_stop - ex.x0)/ex.dx)
    j_0 = int((y_start - ex.y0)/ex.dy)
    j_f = int((y_stop - ex.y0)/ex.dy)
    return grid_x[i_0:i_f], grid_y[j_0:j_f]
        
