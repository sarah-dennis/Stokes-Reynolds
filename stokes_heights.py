# -*- coding: utf-8 -*-
"""
Created on Tue May 21 15:15:15 2024

@author: sarah
"""
import numpy as np
import graphics
from domain import Space

class triangle(Space):
    def __init__(self, x0, xf, y0, yf, N, U, Q, Re, filestr, slopes, wavelen):

        super().__init__(x0, xf, y0, yf, N, U, Q, Re, filestr)
        self.slope_A = slopes[0]
        self.slope_B = slopes[1]
        self.apex = self.Nx//(wavelen*2)    #prev // = 2
        self.spacestr = "Triangle-cavity $Re=%.5f$"%Re
        self.set_space(self.make_space())
        
    def make_space(self):
        grid = np.zeros((self.Ny,self.Nx))
        
        for i in range(self.Nx):
            for j in range(self.Ny):
                # boundary : 0
                # interior : 1
                # exterior : -1
                
                if j == self.Ny-1: # top boundary (lid-driven)
                    grid[j,i] = 0
                    
                elif i < self.apex: # 
                    hj = self.slope_A*i + (self.Ny-1) # yj = h(xi)

                    if j < hj:
                        grid[j,i] = -1
                    elif j > hj:
                        grid[j,i] = 1
                    else:
                        grid[j,i] = 0
                
                elif i > self.apex:
                    hj = self.slope_B * (i - self.apex) # yj = h(xi)
    
                    if j < hj:
                        grid[j,i] = -1
                    elif j > hj:
                        grid[j,i] = 1
                    else:
                        grid[j,i] = 0
                elif i == self.apex and j > 0:
                    grid[j,i] = 1
                
                else:
                    grid[j,i] = 0
        graphics.plot_contour_mesh(grid, self.xs, self.ys, 'space',['space', 'x', 'y'])

        return grid
    
    def stream_interp_EW(self, j, i, psi_k):
        if i < self.apex:
            scale = 1 - (j % self.slope_A)/self.slope_A
        else:
            scale = 1 - (j%self.slope_B)/self.slope_B
        
        return -scale * psi_k

    def stream_interp_S(self, j, i, psi_k):
        if i < self.apex:
            scale = 1 - (i % (1/self.slope_A))*self.slope_A
        else:
            scale = 1 - (i % (1/self.slope_B))*self.slope_B
        return -scale * psi_k

    def streamInlet(self, j):
        if j == self.Ny-1:
            return self.flux
        else:
            return 0
    
    def streamOutlet(self, j):
        if j == self.Ny-1:
            return self.flux
        else:
            return 0
        
    
    def velInlet(self,j):
        if j == self.Ny-1:
            return self.U
        else:
            return 0
        
    def velOutlet(self,j):
        if j == self.Ny-1:
            return self.U
        else:
            return 0
#------------------------------------------------------------------------------

class step(Space):
    def __init__(self, x0, xf, y0, yf, N, U, Q, Re, filestr, x_step, y_step):
        super().__init__(x0, xf, y0, yf, N, U, Q, Re, filestr)
        
        self.i_step = self.Nx//x_step
        self.jf_in = self.Ny//y_step 
        self.jf_out = 0
        self.spacestr = "Backward Facing Step $Re=%.5f$"%Re
        self.set_space(self.make_space())
        
        # constants BCs on velocity, stream, flux 
        self.hf_in = self.y0 + self.jf_in*self.dy
        self.hf_out = self.y0 + self.jf_out*self.dy
        self.H_in = self.yf - self.hf_in
        self.H_out = self.yf - self.hf_out
        
        self.dp_in =  (self.flux - 0.5*self.U*self.H_in) * (-12 / self.H_in**3)
        self.dp_out = (self.flux - 0.5*self.U*self.H_out) * (-12 / self.H_out**3)
        
    def make_space(self):
        grid = np.zeros((self.Ny, self.Nx))
        for j in range(self.Ny):
            for i in range(self.Nx):
                
                if j == self.Ny-1 or i==self.Nx-1:
                    grid[j,i] = 0
                    
                elif j > self.jf_in:
                    if i == 0:
                        grid[j,i] = 0
                    else:
                        grid[j,i] = 1
                elif j == self.jf_in:
                    if i <= self.i_step:
                        grid[j,i] = 0
                    else:
                        grid[j,i] = 1
                    
                else:
                    if i < self.i_step:
                        grid[j,i] = -1
                    elif i == self.i_step or j == 0:
                        grid[j,i] = 0
                    else:
                        grid[j,i] = 1
                        
        return grid
                

    # Boundary conditions on stream and velocity
    def streamInlet(self, j):
        if j > self.jf_in:
            y = self.y0 + j*self.dy
            u_term = self.U* (0.5*(y**2 - self.yf**2) - self.hf_in*(y-self.yf))/self.H_in
            dp_term = -0.5*self.dp_in*( (-1/3)*(y**3 -self.yf**3) + 0.5*(self.yf+self.hf_in)*(y**2-self.yf**2) - self.yf*self.hf_in*(y-self.yf))
            psi = u_term + dp_term + self.flux 
            return psi
        else:
            return 0
    
    def streamOutlet(self, j):
        if j > self.jf_out:
            y = self.y0 + j*self.dy
            u_term = self.U* (0.5*(y**2 - self.yf**2) - self.hf_out*(y-self.yf))/self.H_out
            dp_term = -0.5*self.dp_out*( (-1/3)*(y**3 -self.yf**3) + 0.5*(self.yf+self.hf_out)*(y**2-self.yf**2) - self.yf*self.hf_out*(y-self.yf))
            return u_term + dp_term + self.flux
        else:
            return 0
        
    
    def velInlet(self,j):
        if j > self.jf_in:
            y = self.y0 + j*self.dy

            u = (self.U/self.H_in - 0.5*self.dp_in*(self.yf-y)) * (y-self.hf_in) 

            return  u
        else: 
            return 0
        
    def velOutlet(self,j):
        if j > self.jf_out:
            y = self.y0 + j*self.dy
            u = (self.U/self.H_out - 0.5*self.dp_out*(self.yf-y)) * (y-self.hf_out) 
            return u
        else: 
            return 0
        
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
        self.jf_in = (y_peaks[0][1]-y0)/ self.dy # ys[jf_in] : inlet j-min
        self.jf_out =  (y_peaks[-1][0]-y0)/ self.dy # ys[jf_out] : outlet j-min
        
        self.spacestr = "Backward Facing Step $Re=%.5f$"%Re
        self.set_space(self.make_space())
                                      
        # constants BCs on velocity, stream, flux 
        self.hf_in = y_peaks[0] # hf < h0 measured from y0
        self.hf_out = y_peaks[-1]
        self.H_in = yf - y_peaks[0][1]
        self.H_out = yf - y_peaks[-1][0]
        
        self.dp_in =  (self.flux - 0.5*self.U*self.H_in) * (-12 / self.H_in**3)
        self.dp_out = (self.flux - 0.5*self.U*self.H_out) * (-12 / self.H_out**3)
    
    # 0 : boundary, -1: exterior, 1: interior
    def make_space(self):
        N_regions = len(self.x_peaks)-1
        self.slope = np.zeros(N_regions)
        
        for k in range(N_regions):
            dy = self.y_peaks[k+1][0] - self.y_peaks[k][1]
            dx = self.x_peaks[k+1] - self.x_peaks[k]
            self.slope[k] = dy/dx
        # print(self.slope)
        
        
        grid = np.zeros((self.Ny, self.Nx))
        
        reg = 0 #pwl region
        i_ref = 0 #left index of region 
        
        for i in range(self.Nx):
            
            if self.xs[i] > self.x_peaks[reg+1]:
                reg += 1
                i_ref = i
                
                
            hj = self.slope[reg]*(i - i_ref)*self.dx + self.y_peaks[reg][1] # yj = h(xi)
            # print(hj)
            for j in range(self.Ny):
                if j == self.Ny-1: #upper boundary
                    grid[j,i] = 0
                    
                elif i == 0 and j >= self.jf_in: #inlet boundary
                    grid[j,i] = 0
                    
                elif i == self.Nx-1 and j>= self.jf_out: # outlet boundary
                    grid[j,i] = 0
                else:
                    yj = self.ys[j]
                    if yj < hj:
                        grid[j,i] = -1
                    elif yj > hj:
                        grid[j,i] = 1
                    else:
                        grid[j,i] = 0
        
        
        graphics.plot_contour_mesh(grid, self.xs, self.ys, 'space',['space', 'x', 'y'])
        return grid



    # Boundary conditions on stream and velocity
    def streamInlet(self, j):
        if j > self.jf_in:
            y = self.y0 + j*self.dy
            u_term = self.U* (0.5*(y**2 - self.yf**2) - self.hf_in*(y-self.yf))/self.H_in
            dp_term = -0.5*self.dp_in*( (-1/3)*(y**3 -self.yf**3) + 0.5*(self.yf+self.hf_in)*(y**2-self.yf**2) - self.yf*self.hf_in*(y-self.yf))
            psi = u_term + dp_term + self.flux 
            return psi
        else:
            return 0
    
    def streamOutlet(self, j):
        if j > self.jf_out:
            y = self.y0 + j*self.dy
            u_term = self.U* (0.5*(y**2 - self.yf**2) - self.hf_out*(y-self.yf))/self.H_out
            dp_term = -0.5*self.dp_out*( (-1/3)*(y**3 -self.yf**3) + 0.5*(self.yf+self.hf_out)*(y**2-self.yf**2) - self.yf*self.hf_out*(y-self.yf))
            return u_term + dp_term + self.flux
        else:
            return 0
        
    
    def velInlet(self,j):
        if j > self.jf_in:
            y = self.y0 + j*self.dy

            u = (self.U/self.H_in - 0.5*self.dp_in*(self.yf-y)) * (y-self.hf_in) 

            return  u
        else: 
            return 0
        
    def velOutlet(self,j):
        if j > self.jf_out:
            y = self.y0 + j*self.dy
            u = (self.U/self.H_out - 0.5*self.dp_out*(self.yf-y)) * (y-self.hf_out) 
            return u
        else: 
            return 0
        
        #TODO: make this region interpreting based on i -- save the i-peaks in preprocess
    # def stream_interp_EW(self, j, i, psi_k):
    #     if i < self.apex:
    #         scale = 1 - (j % self.slope_A)/self.slope_A
    #     else:
    #         scale = 1 - (j%self.slope_B)/self.slope_B
        
    #     return -scale * psi_k

    # def stream_interp_S(self, j, i, psi_k):
    #     if i < self.apex:
    #         scale = 1 - (i % (1/self.slope_A))*self.slope_A
    #     else:
    #         scale = 1 - (i % (1/self.slope_B))*self.slope_B
    #     return -scale * psi_k