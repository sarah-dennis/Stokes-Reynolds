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

colour_map_stream = 'plasma'
my_cmap = pp.cm.plasma(np.arange(pp.cm.plasma.N))
my_cmap[:,0:3] *= 1
colour_map_stream = colors.ListedColormap(my_cmap)


# colour_map_mesh = 'plasma'
# colour_map_mesh='PiYG'
# colour_map_mesh='RdYlBu_r'
# colour_map_mesh = 'Spectral_r'
my_cmap = pp.cm.RdYlBu_r(np.arange(pp.cm.RdYlBu_r.N))
# my_cmap = pp.cm.PRGn(np.arange(pp.cm.PRGn.N))
# my_cmap = pp.cm.Spectral_r(np.arange(pp.cm.Spectral_r.N))
my_cmap[:,0:3] *= 0.95
colour_map_mesh = colors.ListedColormap(my_cmap)

# colour_bar_scale=0.015
# colour_bar_scale=0.024
colour_bar_scale=0.04


contour_width = 0.25
stream_width=1

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


#------------------------------------------------------------------------------
def plot_2D(fs, xs, title, axis_labels, color='darkmagenta'):
    fig = pp.figure()
    pp.rcParams['figure.dpi'] = 300
    # marker = None
    marker = 'o'
    pp.plot(xs, fs, color=color, linewidth=.8, marker=marker)

    pp.title(title, fontweight="bold")
    
    pp.xlabel(axis_labels[0])
    
    pp.ylabel(axis_labels[1])
    # pp.ylim(0, 1.1*max(fs))
    pp.minorticks_on()
    return fig

def scatter_2D(fs, xs, title, axis):
    fig = pp.figure()
    pp.scatter(xs, fs, color='b')

    pp.title(title)
    
    pp.xlabel(axis[0])
    
    pp.ylabel(axis[1])
    pp.ylim(0, 1.1*max(fs))
    
    return fig

def plot_2D_inf(fs, xs, title, ax_labels):
    pp.rcParams['figure.dpi'] = 300

    fig, (ax1,ax2) = pp.subplots(1,2,sharey=True, width_ratios=[len(fs)-1,1])
    fig.subplots_adjust(wspace=0.1)
    
    ax1.plot(xs, fs, color='b', linewidth=.8, marker='o')
    ax2.scatter(xs, fs, color='b', linewidth=.8, marker='o')
    
    
    buffer=(xs[1]-xs[0])
    ax1.set_xlim(max(0,xs[0]-buffer),xs[-2]+buffer)
    ax2.set_xlim(xs[-1]-buffer,xs[-1]+buffer)

    
    ax1.spines.right.set_visible(False)
    ax2.spines.left.set_visible(False)
    ax2.yaxis.set_ticks_position('none')
    
    ax2.set_xticks([xs[-1]],['$\infty$\n BFS'])

    pp.sca(ax1)
    pp.title(title, fontweight="bold")
    pp.xlabel(ax_labels[0])
    
    pp.ylabel(ax_labels[1])

    return fig


def plot_2D_inf_multi(fs, xs, title, ax_labels, f_labels ):
    pp.rcParams['figure.dpi'] = 300

    fig, (ax1,ax2) = pp.subplots(1,2,sharey=True, width_ratios=[len(fs[0])-1,1])
    fig.subplots_adjust(wspace=0.1)
    colors = ['forestgreen', 'darkmagenta', 'darkorgange']
    for k in range(len(fs)):
        ax1.plot(xs, fs[k], color=colors[k], linewidth=.8, marker='o',label=f_labels[k])
        ax2.scatter(xs, fs[k], color=colors[k])
    
    
    buffer=(xs[1]-xs[0])
    ax1.set_xlim(max(0,xs[0]-buffer),xs[-2]+buffer)
    ax2.set_xlim(xs[-1]-buffer,xs[-1]+buffer)

    
    ax1.spines.right.set_visible(False)
    ax2.spines.left.set_visible(False)
    ax2.yaxis.set_ticks_position('none')
    
    ax2.set_xticks([xs[-1]],['$\infty$\n BFS'])

    pp.sca(ax1)
    pp.title(title, fontweight="bold")
    pp.xlabel(ax_labels[0])
    
    pp.ylabel(ax_labels[1])
    fig.legend()
    return fig


def plot_2D_multi(fs, xs, title, fun_labels, ax_labels, loc='upper'):
    fig = pp.figure()
    pp.rcParams['figure.dpi'] = 1000
    ax = fig.add_subplot()
    colors = ['r','b','forestgreen', 'darkmagenta', 'darkorgange']
    markers = ['D', 'o', 's', '*', 'H', 'X']
    for i in range(len(fs)):
        
        ax.plot(xs, fs[i], label=fun_labels[i], color=colors[i], linewidth=0.8, marker=markers[i])
    
    #ax.set_xlim([0, 1])
    #ax.set_ylim([0, 1])

    ax.set_xlabel(ax_labels[0])
    ax.set_ylabel(ax_labels[1])
    pp.title(title,  fontweight ="bold")
    pp.minorticks_on()
    if loc== 'upper':
        
        fig.legend(bbox_to_anchor=(0.9, 0.875))
    else:
        fig.legend(bbox_to_anchor=(0.9, 0.275))
    return fig

def plot_2D_multi_multi(fs, xs, title, fun_labels, ax_labels, loc):
    fig = pp.figure()
    pp.rcParams['figure.dpi'] = 1000
    ax = fig.add_subplot()
    colors = ['forestgreen', 'darkmagenta', 'orange', 'firebrick', 'royalblue', 'crimson']
    
    colors = ['firebrick', 'mediumblue','crimson','royalblue','indianred',  'dodgerblue']
    markers = ['D', 'o', 's', '^', 'd', 'h']
    for i in range(len(fs)):
        
        ax.plot(xs[i], fs[i], label=fun_labels[i], color=colors[i], linewidth=0.8, marker=markers[i])
    
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

    stream_density=[1,1]
    pp.streamplot(xs, ys, vx, vy, stream_density, linewidth=0.5, color='k', broken_streamlines=False)
    
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
    ax.set_ylim(0,1.01*max(ys))
    pp.show()
    
        
def plot_stream_heat(vx, vy, xs, ys, color_map, title, ax_labels, vmin, vmax, log_cmap=False, linthresh=1e-18):
    
    pp.rcParams['figure.dpi'] = 1000
    
    pp.figure()
    
    X, Y = np.meshgrid(xs, ys)
    
    stream_density=[1,1]
    
    if log_cmap:
        norm_symLog = colors.AsinhNorm(linthresh, vmin=vmin, vmax=vmax, clip=False)
        stream_plot=pp.streamplot(xs, ys, vx, vy, stream_density, broken_streamlines=False, linewidth=stream_width, color=color_map, cmap=colour_map_stream, norm=norm_symLog)
    else:
        no_norm = colors.CenteredNorm(vcenter=vmin + vmax/2, halfrange=vmax/2)
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
    ax.set_ylim(0,max(ys))
    pp.minorticks_on()
    pp.show()

#------------------------------------------------------------------------------       

def plot_contour(zs, xs, ys, title, labels, log_cmap=False, linthresh=1e-16):
    pp.rcParams["lines.linewidth"] = .5
    pp.rcParams['figure.dpi'] = 1000

    pp.figure()
    
    X, Y = np.meshgrid(xs, ys)
    n_contours = max(zs.shape)//4

    if log_cmap:
        norm_symLog = colors.AsinhNorm(linthresh)#, vmin=-1, vmax=1, clip=True)
        contour_plot = pp.contour(X, Y, zs,  n_contours, cmap='plasma', norm=norm_symLog)
    else:
        contour_plot = pp.contour(X, Y, zs,  n_contours, cmap='plasma')
        
        
    pp.title(title, fontweight="bold")
    pp.xlabel(labels[1])
    pp.ylabel(labels[2])
    pp.colorbar(contour_plot, label=labels[0])
    
    # ax = pp.gca()
    # ax.set_aspect('equal')
    pp.show()



def plot_contour_mesh(zs, xs, ys, title, labels, vmin, vmax, log_cmap=False, linthresh=1e-16, n_contours=20):
    pp.rcParams['figure.dpi'] = 1000
    pp.figure()
    X, Y = np.meshgrid(xs, ys)

    if log_cmap:
        norm_symLog = colors.AsinhNorm(linthresh, vmin=vmin, vmax=vmax, clip=False)
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
    ax.set_ylim(min(ys),max(ys))
    # ax.set_facecolor('black')
    pp.show()    
#------------------------------------------------------------------------------   

def plot_log(fs, xs, title, ax_labels):
    fig = pp.figure()
    pp.rcParams['figure.dpi'] = 300
    pp.loglog(xs, fs, color='b')   

    pp.title(title)
    
    pp.xlabel(ax_labels[0])
    pp.ylabel(ax_labels[1])
    
    return fig

def plot_log_x_multi(fs, xs, title, f_labels,ax_labels):
    fig = pp.figure()
    pp.rcParams['figure.dpi'] = 300
    pp.xscale('log')
    colors = ['firebrick', 'royalblue', 'forestgreen', 'darkorchid', 'goldenrod']
    ax = fig.gca()
    for i in range(len(fs)):
        ax.plot(xs, fs[i], label=f_labels[i], color=colors[i], marker='x', markevery=1)
    
    pp.title(title)
    
    pp.xlabel(ax_labels[0])
    pp.ylabel(ax_labels[1])
    
    pp.title(title,  fontweight ="bold")
    fig.legend()
    
    
    return fig

    
def plot_log_multi(fs, xs, title, f_labels, ax_labels, linthresh=1e-16, O1=1e-2,O2=1e-3):
    pp.rcParams['figure.dpi'] = 300
    fig = pp.figure()
    
    
    ax = fig.add_subplot()
    colors = [ 'forestgreen', 'darkorchid', 'darkorange', 'royalblue', 'firebrick']    
    markers = ['D', 'o', 's', '*', 'H', 'X']

    pp.rcParams["lines.linewidth"] =1
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
    fig.legend(bbox_to_anchor=(0.3, 0.425))
    
    return fig

