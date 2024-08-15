# -*- coding: utf-8 -*-
"""
Created on Tue May 21 15:15:15 2024

@author: sarah
"""
import numpy as np
from domain import Space

class triangle(Space):
    def __init__(self, x0, xf, y0, yf, N, U, Q, Re, filestr):

        super().__init__(x0, xf, y0, yf, N, U, Q, Re, filestr)
        # self.slope = int(2*yf/xf)
        self.slope_A = -4
        self.slope_B = 4
        self.apex = self.Nx//2
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