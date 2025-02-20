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

            y=height.ys[j]
            
            i = 0 #inlet bc
            if y <= height.hs[i]:
                ps_adj[j,i]=reyn_ps[i]
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
                ps_adj[j,i]=reyn_ps[i]
            else:
                ps_adj[j,i]=None 
                
            i=height.Nx-2 #outlet bkwd difference
            if y <= height.hs[i]:   
               
                dP_E = (3*reyn_ps[i+1]-4*reyn_ps[i]+reyn_ps[i-1])/(2*height.dx) 
                dP_W = (reyn_ps[i]-reyn_ps[i-2])/(2*height.dx)
                
                adj_E=(hs[i+1]-y)*dP_E/2 + height.U*height.visc/hs[i+1]
                adj_W=(hs[i-1]-y)*dP_W/2 + height.U*height.visc/hs[i-1]

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
                    h_EE = height.hs[i+2] 
                    h_WW = height.hs[i-2]
                    
                    if y <= h_EE and y <= h_WW: # interior 

                        dP_E = (reyn_ps[i+2]-reyn_ps[i])/(2*height.dx)
                        dP_W = (reyn_ps[i]-reyn_ps[i-2])/(2*height.dx)
                    
                    elif y > h_EE and y <= h_WW: # east oob 
                        
                        dP_W = (reyn_ps[i]-reyn_ps[i-2])/(2*height.dx)

                        if y <= h_E: #
                            dP_E = (reyn_ps[i+1]-reyn_ps[i])/height.dx 
                        else:
                            dP_E = dP_W
                            
                    elif y <= h_EE and y > h_WW: #west oob
                        
                        dP_E = (reyn_ps[i+2]-reyn_ps[i])/(2*height.dx) 
                        
                        if y <= h_W:
                            dP_W = (reyn_ps[i] - reyn_ps[i-1])/height.dx 
                        else:
                            dP_W = dP_E
                            
                    else: #both east and west nbr oob 
                        dP_E = 0
                        dP_W = 0
                    
                    adj_E = (h_E-y) * dP_E/2 - height.U*height.visc/h_E
                    adj_W = (h_W-y) * dP_W/2 - height.U*height.visc/h_W

                    adj = y*(adj_E - adj_W)/(2*height.dx)

                    ps_adj[j,i] = reyn_ps[i] + adj
    
        # ps_adj = np.flip(ps_adj, axis=0)

        return ps_adj
    

