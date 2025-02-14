 #!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec  5 13:28:42 2024

@author: sarahdennis
"""
import numpy as np

class Pressure:
    def __init__(self, height, ps):
        self.ps_1D = ps
        self.ps_2D = self.make_2D_ps(height,ps)
        
    def make_2D_ps(s_elf,height,ps):
        ps_2D = np.zeros((height.Ny, height.Nx))
        
        for i in range(height.Nx):
            for j in range(height.Ny):
                
                y = height.ys[j]
                if y <= height.hs[i]:
                    ps_2D[j,i] = ps[i]
                else:
                    ps_2D[j,i] = None
                    
                
        ps_2D = np.flip(ps_2D, axis=0)
        return ps_2D
            
    
class Adj_Pressure:
    
    def __init__(self, height, reyn_ps):
        self.ps_2D= self.make_adj_ps(height, reyn_ps)
        
    def make_adj_ps(self,height,reyn_ps):
        ps_adj = np.zeros((height.Ny, height.Nx))
        hs = height.hs
        
    
        for j in range(height.Ny):
            ps_adj[j,0]=reyn_ps[0]
            ps_adj[j,height.Nx-1]=reyn_ps[height.Nx-1]
            
            y=height.ys[j]
            
            i=1                    
            dP_E = (reyn_ps[i+2]-reyn_ps[i])/(2*height.dx)
            dP_W = (-3*reyn_ps[i-1]+4*reyn_ps[i]-reyn_ps[i+1])/(2*height.dx) 
            
            adj_E=(hs[i+1]-y)*dP_E/2 + height.U*height.visc/hs[i+1]
            adj_W=(hs[i-1]-y)*dP_W/2 + height.U*height.visc/hs[i-1]
            adj = y*(adj_E - adj_W)/(2*height.dx)
            ps_adj[j,i] = reyn_ps[i] + adj

            i=height.Nx-2
            dP_E = (3*reyn_ps[i+1]-4*reyn_ps[i]+reyn_ps[i-1])/(2*height.dx) 
            dP_W = (reyn_ps[i]-reyn_ps[i-2])/(2*height.dx)
            
            adj_E=(hs[i+1]-y)*dP_E/2 + height.U*height.visc/hs[i+1]
            adj_W=(hs[i-1]-y)*dP_W/2 + height.U*height.visc/hs[i-1]

            adj = y*(adj_E - adj_W)/(2*height.dx)

            ps_adj[j,i] = reyn_ps[i] + adj
            
        
        for i in range(2,height.Nx-2):
            for j in range(height.Ny):
                
                y = height.ys[j]
                if y > hs[i]:
                    ps_adj[j,i] = None
                else:
                    
                    dP_E = (reyn_ps[i+2]-reyn_ps[i])/(2*height.dx)
                    dP_W = (reyn_ps[i]-reyn_ps[i-2])/(2*height.dx)
                    
                    adj_E=(hs[i+1]-y)*dP_E/2 + height.U*height.visc/hs[i+1]
                    adj_W=(hs[i-1]-y)*dP_W/2 + height.U*height.visc/hs[i-1]

                    adj = y*(adj_E - adj_W)/(2*height.dx)

                    ps_adj[j,i] = reyn_ps[i] + adj
    
        ps_adj = np.flip(ps_adj, axis=0)

        return ps_adj
    

