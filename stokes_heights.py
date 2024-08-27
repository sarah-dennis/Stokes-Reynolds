# -*- coding: utf-8 -*-
"""
Created on Tue May 21 15:15:15 2024

@author: sarah
"""
import numpy as np
import graphics
import math
from domain import Space

#------------------------------------------------------------------------------
#TODO
class slider(Space):
    # peak_xs = [x0, ..., xi,..., xf] : x0 < xi < xf
    # peak_ys = [(0,h_in),...,(hi_left, hi_right),...,(h_out,0)]
    
    
    def __init__(self, x0, xf, y0, yf, N, U, Q, Re, filestr, x_peaks, y_peaks):
        super().__init__(x0, xf, y0, yf, N, U, Q, Re, filestr)
        
        # self.i_step = self.Nx// (peak_xs[1]-x0)
        self.x_peaks=x_peaks
        self.y_peaks = y_peaks
        
        #peaks must have integer indices in the grid! 
        # self.jf_in = y_peaks[0][1]/ self.dy # ys[jf_in] : inlet j-min
        # self.jf_out =  y_peaks[-1][0] / self.dy # ys[jf_out] : outlet j-min
        
        self.spacestr = "Textured Slider $Re=%.4f$"%Re
        self.set_space(self.make_space())
                                      
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
        
        self.N_regions = len(self.x_peaks)-1
        self.slopes = np.zeros(self.N_regions)
        for k in range(self.N_regions):
            dy = self.y_peaks[k+1][0] - self.y_peaks[k][1]
            dx = self.x_peaks[k+1] - self.x_peaks[k]
            self.slopes[k] = dy/dx
        
        grid = np.zeros((self.Ny, self.Nx))
    
        reg = 0 #pwl region        
        i_ref = 0
        for i in range(self.Nx):
            if self.xs[i] == self.x_peaks[reg]:
                h_left = self.y_peaks[reg][0]
                h_right = self.y_peaks[reg][1] 
                i_ref = i
                reg +=1
                
            else:
                h = self.slopes[reg-1]*(i - i_ref)*self.dx + self.y_peaks[reg-1][1]

            for j in range(self.Ny):
                y = self.ys[j]
                if j == self.Ny-1: #upper boundary
                    grid[j,i] = 0
                    
                elif i == 0: #inlet boundary
                    if y >= self.y_peaks[0][1]: 
                        grid[j,i] = 0
                    else:
                        grid[j,i]=-1
                    
                elif i == self.Nx-1:# outlet boundary
                    if y >= self.y_peaks[-1][0]: 
                        grid[j,i] = 0
                    else:
                        grid[j,i]=-1
                
                elif i == i_ref: # pwl region change              
                    
                    if y == h_left or y == h_right: # true boundary point
                        grid[j,i] = 0
                        
                    elif h_left < h_right:
                        if h_left < y and y < h_right: # x=h(y) boundary
                            grid[j,i] = 0
                        elif y > h_right:
                            grid[j,i]=1
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
        return grid



    # Boundary conditions on stream and velocity
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
      
# Interior x_ij = x_k has exterior nbr x_st

    def interp_E_W(self, t, s, psi_k):
        x = self.xs[s]
        for k in range(self.N_regions):
            if x >= self.x_peaks[k]:
                slope = self.slopes[k]

        if slope == 0:
            scale = 0
        else:
            scale = 1 - (t%slope)/slope
        
        return -scale * psi_k

    def interp_S(self, t, s, psi_k):
        x = self.xs[s]
        for k in range(self.N_regions):
            if x >= self.x_peaks[k]:
                slope = self.slopes[k]

        scale = 1 - (s % (1/slope))*slope
        return -scale * psi_k

    def interp_NE_NW(self, t, s, psi_N):
        x = self.xs[s]
        for k in range(self.N_regions):
            if x >= self.x_peaks[k]:
                slope = self.slopes[k]

        scale = 1 - (t%slope)/slope
        
        return -scale * psi_N

    def interp_SE_SW(self, t, s, psi_EW, psi_S):
        x = self.xs[s]
        for k in range(self.N_regions):
            if x >= self.x_peaks[k]:
                slope = self.slopes[k]

        if slope == 0:
            return 0
        
        elif abs(slope) < 1:    
            scale = 1 - (s%(1/slope))*slope
            return -scale * psi_EW
        
        else:
            scale = 1 - (t % slope)/slope
            return -scale * psi_S

    

