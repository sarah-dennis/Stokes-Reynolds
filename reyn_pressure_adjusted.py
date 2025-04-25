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
    dx= height.dx
    for j in range(height.Ny):

        y = height.ys[j]
        
        i = 0 #inlet bc
        if y <= height.hs[i]:

            dP_E = dm.center_first(dx, reyn_ps[i : i+3]) # centered at i+1 = 1
            dP_W = dP_E# 0 # i-1 < 0

            adj_E=(hs[i+1]-y)*dP_E/2 + height.U*height.visc/hs[i+1] 
            adj_W=(hs[i]-y)*dP_W/2 + height.U*height.visc/hs[i]
            adj = y*(adj_E - adj_W)/(2*height.dx)
            
            ps_adj[j,i] = reyn_ps[i] + adj
        else:
            ps_adj[j,i]=None  
            
        i=1 # inlet fwd difference   
        if y <= height.hs[i]:   
            
            dP_E = dm.center_first(dx, reyn_ps[i : i+3]) # centered at i+1 = 2
            
            dP_W = dm.right_first(dx, reyn_ps[i-1 : i+2]) # right from i-1 = 0

            adj_E=(hs[i+1]-y)*dP_E/2 + height.U*height.visc/hs[i+1] 
            adj_W=(hs[i-1]-y)*dP_W/2 + height.U*height.visc/hs[i-1]
            adj = y*(adj_E - adj_W)/(2*height.dx)
            ps_adj[j,i] = reyn_ps[i] + adj
            
        else:
            ps_adj[j,i] = None
            
        i = height.Nx-1 # outlet bc
        if y <= height.hs[i]:
            
            
            dP_W = dm.center_first(dx, reyn_ps[i-2 : i+1]) # centered at i-1 = Nx-2
            dP_E = dP_W # i + 1 > Nx-1
            
            adj_E = (hs[i]-y)*dP_E/2 + height.U*height.visc/hs[i]
            adj_W = (hs[i-1]-y)*dP_W/2 + height.U*height.visc/hs[i-1]

            adj = y*(adj_E - adj_W)/(2*height.dx)
            ps_adj[j,i] = reyn_ps[i] + adj
        else:
            ps_adj[j,i]=None 
            
        i=height.Nx-2 #outlet bkwd difference
        if y <= height.hs[i]:   
           
            dP_E = dm.left_first(dx, reyn_ps[i-1 : i+2]) #left from i+1 = Nx-1
            dP_W = dm.center_first(dx, reyn_ps[i-2 : i+1]) #centered at i-1 = Nx-3
            
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
                   
                dP_E = dm.center_first(dx, reyn_ps[i : i+3])
                dP_W = dm.center_first(dx, reyn_ps[i-2 : i+1])
                
                adj_E = (h_E-y) * dP_E/2 - height.U*height.visc/h_E
                adj_W = (h_W-y) * dP_W/2 - height.U*height.visc/h_W

                adj = y*(adj_E - adj_W)/(2*height.dx)

                ps_adj[j,i] = reyn_ps[i] + adj

    # ps_adj = np.flip(ps_adj, axis=0)
    # graphics.plot_2D(ps_adj[0,:], height.xs, 'p(x,y0)', ['p', 'x', 'y'])
    return ps_adj 