# -*- coding: utf-8 -*-
"""
Created on Tue May 21 15:15:15 2024

@author: sarah
"""
import numpy as np
from domain import Height

class triangle(Height):
    def __init__(self, x0, xf, y0, yf, N, U, Re, filestr):
        # N even --> (xL-x0)/2 is triangle vertex
        # slope dividing 2N  --> maximal true boundary points

        self.slope = int(2*yf/xf)
        self.filename = filestr

        Nx = (xf-x0)*N + 1

        self.apex = Nx//2
        hs = np.zeros(Nx) # Reynolds needs this, Stokes its built into dPsi
        # TODO: make hs since it could be helpful to make for plotting

        super().__init__(x0, xf, y0, yf, N, hs, U, Re, filestr)
        self.bndry_nbrs = self.boundary_helper()
        
    def is_interior(self, i, j):
        if i == 0 or j == 0 or i == self.Nx-1 or j == self.Ny-1:
            return False 
   
        elif  (i < self.apex and j <= (self.Ny-1) - self.slope*i):
            return False
        
        elif  (i > self.apex and j <= self.slope * (i - self.apex)):
            return False
          
        else:
            return True
 
# used in update_rhs        
    def boundary_helper(self):
        bndry_nbrs= np.zeros((self.Ny, 6))
        
        for j in range(self.slope, self.Ny): 
            dj = j % self.slope
            di = int(j//self.slope) 
            
            if dj != 0 : # row j has no true boundary values
            
                # i-index of interior points in row j who may have external nbrs
                i_left = self.apex - di
                i_right = self.apex + di
                
                # North West
                bndry_nbrs[j][0] = not self.is_interior(i_left-1, j+1) 
                # West
                bndry_nbrs[j][1] = not self.is_interior(i_left-1, j)
                # South West
                bndry_nbrs[j][2] = not self.is_interior(i_left-1, j-1)  
    

                # South East
                bndry_nbrs[j][3] = not self.is_interior(i_right+1, j-1)
                
                #East
                bndry_nbrs[j][4] = not self.is_interior(i_right+1, j)

                #North East
                bndry_nbrs[j][5] = not self.is_interior(i_right+1, j+1)

        return bndry_nbrs

#------------------------------------------------------------------------------

class step(Height):
    def __init__(self, x0, xf, y0, yf, N, U, Q, Re, filestr, x_step, y_step):
        #y0=hmin, yf=hmax
        #x0=inlet, xf=outlet
                    
        Nx = (xf-x0)*N + 1
        
        Ny = (yf - y0)*N + 1
          
        hs = np.zeros(Nx) #update later
        
        self.filename = filestr 
        super().__init__(x0, xf, y0, yf, N, hs, U, Re, filestr)
        
        #xi where step occurs
        self.i_step = Nx//x_step
    
        #yj = h(xi) lower indices at inlet & outlet
        self.jf_in = Ny//y_step 
        self.jf_out = 0
        
        self.hf_in = self.y0 + self.jf_in*self.dy
        self.hf_out = self.y0 + self.jf_out*self.dy
        self.hs = self.make_hs()
        
        self.H_in = self.yf - self.hf_in
        self.H_out = self.yf - self.hf_out

        self.flux = Q # = stream(x, hc=yf)
        self.dp_in =  (self.flux - 0.5*self.U*self.H_in) * (-12 / self.H_in**3)
        self.dp_out = (self.flux - 0.5*self.U*self.H_out) * (-12 / self.H_out**3)
        
 
        
    def make_hs(self):
        hs = np.zeros(self.Nx)   
        for i in range(self.Nx):
            if i <= self.i_step:
                hs[i] = self.hf_in
            else:
                hs[i] = self.hf_out
        return hs
        
        
#used in update rhs    
    def is_interior(self, i, j):
        if i == 0 or j == 0 or i == self.Nx-1 or j == self.Ny-1:
            return False 
   
        elif  (i <= self.i_step and j <= self.jf_in):
            return False
        
        elif  (i >= self.i_step and j <= self.jf_out):
            return False

        else:
            return True
    
#used in uv approx
    def is_lowerbndry(self, i, j):
        if  (i <= self.i_step and j == self.jf_in):
            return True
        
        elif  (i >= self.i_step and j == self.jf_out):
            return True
        
        elif (i == self.i_step and j < self.jf_in):
            return True
        else:
            return False

    
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