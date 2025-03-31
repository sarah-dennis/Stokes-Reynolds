# -*- coding: utf-8 -*-
"""
Created on Tue Feb 25 15:16:58 2025

@author: sarah
"""
import numpy as np
import domain as dm
import graphics

def make_adj_ps(height, reyn_ps):
    ps_adj = np.zeros((height.Ny, height.Nx))
    hs = height.hs

    for j in range(height.Ny):

        y = height.ys[j]
        
        i = 0 #inlet bc
        if y <= height.hs[i]:
            ps_adj[j,i]=reyn_ps[i]
            
            dP_E = (reyn_ps[i+2]-reyn_ps[i])/(2*height.dx)
            
            dP_W = 0# (reyn_ps[i+1]-reyn_ps[i])/(height.dx) 

            adj_E=(hs[i+1]-y)*dP_E/2 + height.U*height.visc/hs[i+1] 
            adj_W=(hs[i]-y)*dP_W/2 + height.U*height.visc/hs[i]
            adj = y*(adj_E - adj_W)/(2*height.dx)
            ps_adj[j,i] = reyn_ps[i] + adj
        else:
            ps_adj[j,i]=None  
            
        i=1 # inlet fwd difference   
        if y <= height.hs[i]:   
            
            dP_E = (reyn_ps[i+2]-reyn_ps[i])/(2*height.dx)
            
            dP_W = (-3*reyn_ps[i-1]+4*reyn_ps[i]-reyn_ps[i+1])/(2*height.dx) 

            adj_E=(hs[i+1]-y)*dP_E/2 + height.U*height.visc/hs[i+1] 
            adj_W=(hs[i-1]-y)*dP_W/2 + height.U*height.visc/hs[i-1]
            adj = y*(adj_E - adj_W)/(2*height.dx)
            ps_adj[j,i] = reyn_ps[i] + adj
        else:
            ps_adj[j,i] = None
            
        i = height.Nx-1 # outlet bc
        if y <= height.hs[i]:
            
            dP_E = 0 #(reyn_ps[i] -reyn_ps[i-1])/(height.dx) 
            dP_W = (reyn_ps[i] -reyn_ps[i-2])/(2*height.dx)
            
            adj_E = (hs[i]-y)*dP_E/2 + height.U*height.visc/hs[i]
            adj_W = (hs[i-1]-y)*dP_W/2 + height.U*height.visc/hs[i-1]

            adj = y*(adj_E - adj_W)/(2*height.dx)
            ps_adj[j,i] = reyn_ps[i] + adj
        else:
            ps_adj[j,i]=None 
            
        i=height.Nx-2 #outlet bkwd difference
        if y <= height.hs[i]:   
           
            dP_E = (3*reyn_ps[i+1] -4*reyn_ps[i] +reyn_ps[i-1])/(2*height.dx) 
            dP_W = (reyn_ps[i] -reyn_ps[i-2])/(2*height.dx)
            
            adj_E = (hs[i+1]-y)*dP_E/2 + height.U*height.visc/hs[i+1]
            adj_W = (hs[i-1]-y)*dP_W/2 + height.U*height.visc/hs[i-1]

            adj = y*(adj_E - adj_W)/(2*height.dx)
            ps_adj[j,i] = reyn_ps[i] + adj
        else:
            ps_adj[j,i] = None
    
    for i in range(2,height.Nx-2): 
        for j in range(height.Ny):
            
            y = height.ys[j]
            
            if y > hs[i]:     #exterior
                ps_adj[j,i] = None
                
            else:  
                
                #find dP_E = dp/px @ x_i+1 and dP_W = dp/dx # x_i-1
                h_E = height.hs[i+1]
                h_W = height.hs[i-1]
                   
                dP_E = (reyn_ps[i+2]-reyn_ps[i])/(2*height.dx)
                dP_W = (reyn_ps[i]-reyn_ps[i-2])/(2*height.dx)
                
                adj_E = (h_E-y) * dP_E/2 - height.U*height.visc/h_E
                adj_W = (h_W-y) * dP_W/2 - height.U*height.visc/h_W

                adj = y*(adj_E - adj_W)/(2*height.dx)

                ps_adj[j,i] = reyn_ps[i] + adj

    # ps_adj = np.flip(ps_adj, axis=0)
    # graphics.plot_2D(ps_adj[0,:], height.xs, 'p(x,y0)', ['p', 'x', 'y'])
    return ps_adj 