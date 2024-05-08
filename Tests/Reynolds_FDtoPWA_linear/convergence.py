# -*- coding: utf-8 -*-
"""
Created on Tue Mar  7 16:16:13 2023

@author: sarah
"""

import pyplot as pp

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
    
Nxs = [1e2, 1e3, 1e4]

#O(2)
bigO = [1, 1e-2, 1e-4]

#1 linear
errs_1 = [6.264302607166528e-2, 6.267360573772862e-4, 6.257021379241223e-06]

#2 linear
errs_2 = [3.050932554811787e-1, 3.1518906103027433e-3, 3.14755312391668e-05]

#3 linear
errs_3 = [1.0289259062325318e-1,1.0482565594482907e-3,1.0463327988929905e-05]

fs = [bigO, errs_1, errs_2, errs_3]
f_labels = ['$\Theta(N^2)$','linear', '2-linear', '3-linear']
ax_labels = ['$N$ grid points', 'Maximum Error']
title = 'Error in Reynolds Finite Difference to Analytic Pressure \n for Piecewise-Linear Height'
plot_log_multi(fs, Nxs, title, f_labels, ax_labels)


    
