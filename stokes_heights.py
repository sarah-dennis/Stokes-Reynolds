# -*- coding: utf-8 -*-
"""
Created on Tue May 21 15:15:15 2024

@author: sarah
"""
import numpy as np
import math
from domain import Space

#------------------------------------------------------------------------------
class PWLinear(Space):
    # peak_xs = [x0, ..., xi,..., xf] : x0 < xi < xf
    # peak_ys = [(yf,h_in),...,(hi_left, hi_right),...,(h_out,yf)]
    
    def __init__(self, x0, xf, y0, yf, N, U, Q, Re, filestr, x_peaks, y_peaks):
        super().__init__(x0, xf, y0, yf, N, U, Q, Re, filestr)

        self.x_peaks = x_peaks
        self.y_peaks = y_peaks
        self.N_regions = len(self.x_peaks)-1
            
        # peaks must fall on the grid
        self.make_space()
        self.spacestr = "Textured Slider $Re=%.4f$"%Re  
                            
        # constants BCs on velocity, stream, flux 
        self.hf_in = y_peaks[0][1] # hf < h0 measured from y0
        self.hf_out = y_peaks[-1][0]
        self.H_in = yf - y_peaks[0][1]
        self.H_out = yf - y_peaks[-1][0]
        
        if self.H_in == 0 and self.H_out == 0:
            self.dp_in = 0
            self.dp_out = 0 
        else:
            self.dp_in =  (self.flux - 0.5*self.U*self.H_in) * (-12 / self.H_in**3)
            self.dp_out = (self.flux - 0.5*self.U*self.H_out) * (-12 / self.H_out**3)
    
    # 0 : boundary, -1: exterior, 1: interior
    def make_space(self):

        slopes = np.zeros(self.N_regions) #TODO: try to get rid 

        for k in range(self.N_regions):
            dh = self.y_peaks[k+1][0] - self.y_peaks[k][1]
            dx = self.x_peaks[k+1] - self.x_peaks[k]
            slopes[k] = dh/dx
    
        hs = np.zeros(self.Nx)
        
        grid = np.zeros((self.Ny, self.Nx))
    
        reg = 0 #pwl region        
        i_ref = 0
        for i in range(self.Nx):
            if self.xs[i] == self.x_peaks[reg]:
                h_left = self.y_peaks[reg][0]
                h_right = self.y_peaks[reg][1] 
                i_ref = i
                reg +=1
                hs[i] = h_right # choice for 1D
            else:
                h = slopes[reg-1]*(i - i_ref)*self.dx + self.y_peaks[reg-1][1]
                hs[i] = h
            
            for j in range(self.Ny):
                y = self.ys[j]
                if j == self.Ny-1: #upper boundary
                    grid[j,i] = 0
                    
                elif i == 0: #inlet boundary
                    if y >= self.y_peaks[0][1]: 
                        grid[j,i] = 0
                    else:
                        grid[j,i] = -1
                    
                elif i == self.Nx-1:# outlet boundary
                    if y >= self.y_peaks[-1][0]: 
                        grid[j,i] = 0
                    else:
                        grid[j,i] = -1
                
                elif i == i_ref: # pwl region change              
                    
                    if y == h_left or y == h_right: # true boundary point
                        grid[j,i] = 0
                        
                    elif h_left < h_right:
                        if h_left < y and y < h_right: # x=h(y) boundary
                            grid[j,i] = 0
                        elif y > h_right:
                            grid[j,i] = 1
                        else:
                            grid[j,i] = -1
                    
                    else:
                        if h_left > y and y > h_right: # x=h(y) boundary
                            grid[j,i] = 0
                        elif y > h_left:
                            grid[j,i]=1
                        else:
                            grid[j,i] = -1

                else:
                    if math.isclose(y, h): # ~true boundary point ( rational slope*dx)
                        grid[j,i] = 0
                    elif y > h:
                        grid[j,i] = 1
                    else:
                        grid[j,i] = -1
        
        
        # graphics.plot_contour_mesh(grid, self.xs, self.ys, 'space',['space', 'x', 'y'])
        self.space = grid
        self.slopes = slopes
        self.hs = hs

#------------------------------------------------------------------------------
# # Boundary conditions on stream and velocity
#------------------------------------------------------------------------------
   
    def streamInlet(self, j):
        y = self.y0 + j*self.dy
        if y >= self.y_peaks[0][1]:
            y = self.y0 + j*self.dy
            u_term = self.U* (0.5*(y**2 - self.yf**2) - self.hf_in*(y-self.yf))/self.H_in
            dp_term = -0.5*self.dp_in*( (-1/3)*(y**3 -self.yf**3) + 0.5*(self.yf+self.hf_in)*(y**2-self.yf**2) - self.yf*self.hf_in*(y-self.yf))
            psi = u_term + dp_term + self.flux 
            return psi
        else:
            return 0
    
    def streamOutlet(self, j):
        y = self.y0 + j*self.dy
        if y >= self.y_peaks[-1][0]:
            y = self.y0 + j*self.dy
            u_term = self.U* (0.5*(y**2 - self.yf**2) - self.hf_out*(y-self.yf))/self.H_out
            dp_term = -0.5*self.dp_out*( (-1/3)*(y**3 -self.yf**3) + 0.5*(self.yf+self.hf_out)*(y**2-self.yf**2) - self.yf*self.hf_out*(y-self.yf))
            return u_term + dp_term + self.flux
        else:
            return 0
        
    
    def velInlet(self,j):
        y = self.y0 + j*self.dy
        if y >= self.y_peaks[0][1]:
            

            u = (self.U/self.H_in - 0.5*self.dp_in*(self.yf-y)) * (y-self.hf_in) 

            return  u
        else: 
            return 0
        
    def velOutlet(self,j):
        y = self.y0 + j*self.dy
        if y >= self.y_peaks[-1][0]:

            u = (self.U/self.H_out - 0.5*self.dp_out*(self.yf-y)) * (y-self.hf_out) 
            return u
        else: 
            return 0
#------------------------------------------------------------------------------
# Boundary interpolation
#------------------------------------------------------------------------------
 
# x_ij interior node
# x_st exterior nbr node 

    def interp_E_W(self, i,j, s,t, v_ij, v_bdry=0):
        
        slope = (self.hs[i]-self.hs[s])/(self.xs[i]-self.xs[s])
        x_bdry = self.xs[i] + (self.ys[j]-self.hs[i])/slope

        # y_bdry = self.ys[j]
        
        v1 = self.xs[s]-self.xs[i]
        v2 = x_bdry-self.xs[i]
    
        scale =  v2/v1
  
        v_st = v_bdry - v_ij*scale      

        return v_st
    
    def interp_S(self, i,j, s,t, v_ij, v_bdry=0):
        
        #x_bdry = self.xs[i]
        
        y_bdry = self.hs[i] 
        
        v1 = self.ys[t]-self.ys[j] 
        v2 = y_bdry-self.ys[j]
        
        scale = v1/v2 -1

        v_st = v_bdry-v_ij * scale


        return v_st

    def interp_NE(self, i,j, s,t, v_ij, v_bdry=0): 
        
        slope = (self.hs[i]-self.hs[s])/(self.xs[i]-self.xs[s])
        
        x_bdry = (self.ys[j]-self.hs[i])/(slope-1)  +self.xs[i]
        
        y_bdry = (x_bdry-self.xs[i]) +self.ys[j] 
        
        #NE : s > bndry > i, t > bndry > j
        v1 = np.sqrt((self.xs[s]-self.xs[i])**2 + (self.ys[t]-self.ys[j])**2)
        v2 = np.sqrt((x_bdry-self.xs[i])**2 + (y_bdry-self.ys[j])**2)
        
        scale = v1/v2 -1
        
        v_st = v_bdry -v_ij * scale
        return v_st
    
    def interp_SW(self, i,j, s,t, v_ij, v_bdry=0): 
        
        slope = (self.hs[i]-self.hs[s])/(self.xs[i]-self.xs[s])
        
        x_bdry = (self.ys[j]-self.hs[i])/(slope-1)  +self.xs[i]
        
        y_bdry = (x_bdry-self.xs[i]) +self.ys[j] 
        
        #NE : s < bndry < i, t < bndry < j
        v1 = np.sqrt((self.xs[i]-self.xs[s])**2 + (self.ys[j]-self.ys[t])**2)
        v2 = np.sqrt((self.xs[i]-x_bdry)**2 + (self.ys[j]-y_bdry)**2)
        
        scale = v1/v2 -1
        
        v_st = v_bdry-v_ij * scale
        return v_st
    
    
    def interp_NW(self, i,j, s,t, v_ij, v_bdry=0): 

        slope = (self.hs[i]-self.hs[s])/(self.xs[i]-self.xs[s])
        
        x_bdry = (self.ys[j]-self.hs[i])/(slope+1) +self.xs[i]
        y_bdry = -(x_bdry - self.xs[i]) + self.ys[j] 
        
        # NW: s < bndry < i : t > bndry > j 

        v1 = np.sqrt((self.xs[i]-self.xs[s])**2 + (self.ys[t]-self.ys[j])**2)
        v2 = np.sqrt((self.xs[i]-x_bdry)**2 + (y_bdry-self.ys[j])**2)
        
        scale = v1/v2 -1
        
        v_st = v_bdry-v_ij * scale
        return v_st
    
    def interp_SE(self, i,j, s,t, v_ij, v_bdry=0): 

        slope = (self.hs[i]-self.hs[s])/(self.xs[i]-self.xs[s])
        
        x_bdry = (self.ys[j]-self.hs[i])/(slope+1) +self.xs[i]
        y_bdry = -(x_bdry - self.xs[i]) + self.ys[j] 

        #SE : s > bndry > i : t < bndry < j
        v1 = np.sqrt((self.xs[s]-self.xs[i])**2 + (self.ys[j]-self.ys[t])**2)
        v2 = np.sqrt((x_bdry-self.xs[i])**2 + (self.ys[j]-y_bdry)**2)
        
        scale = v1/v2 -1
        
        v_st = v_bdry -v_ij * scale
        return v_st




























