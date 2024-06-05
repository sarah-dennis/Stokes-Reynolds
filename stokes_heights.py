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
        #y0 = 0, yf=2
        self.filename = filestr             
        Nx = (xf-x0)*N + 1
        
        Ny = (yf - y0)*N + 1
          
        hs = np.zeros(Nx)
        
        super().__init__(x0, xf, y0, yf, N, hs, U, Re, filestr)
        

        self.step = Nx//x_step
    
        self.hj_in = Ny//y_step 
        self.hx_in = self.y0 + self.hj_in*self.dy
        
        self.hj_out = 0
        self.hx_out = self.y0 + self.hj_out*self.dy

        self.flux = Q
        
    
    def is_interior(self, i, j):
        if i == 0 or j == 0 or i == self.Nx-1 or j == self.Ny-1:
            return False 
   
        elif  (i <= self.step and j <= self.hj_in):
            return False
        
        elif  (i > self.step and j <= self.hj_out):
            return False

        else:
            return True
        
    def stream(self, y, h0, hx):
        dp = (self.flux - 0.5*self.U*(h0-hx)) * (-12 * (h0-hx)**-3)
        
        v = self.U/(h0-hx) * (0.5*y**2 - hx*y)
        
        p = -0.5*dp*((-1/3)*y**3 + (1/2)*(h0+hx)*y**2 - h0*hx*y)
        
        return v + p# + self.flux
            
    def velocity(self, y, h0, hx):
        dp = (self.flux - (self.U/2)*(h0-hx))*(-12/(h0-hx)**3)
        return self.U/(h0-hx) - 0.5*dp*(h0-y)*(y-hx)
    
    
    def streamInlet(self, y):
        h0 = self.yf           #upper boundary moving y=h0
        hx = self.hx_in        #lower boundary y=h(x0)
        
        if y > hx:

            return self.stream(y, h0, hx)
        else:
            return 0

   
    def streamOutlet(self, y):
        h0 = self.yf            #upper boundary moving y=h0
        hx = self.hx_out        #lower boundary y=h(xf)
        
        if y > hx:
            
            return self.stream(y, h0, hx)
        else:
            return 0
        
        
        
    def velInlet(self,y):
        h0 = self.yf           #upper boundary moving y=h0
        hx = self.hx_in        #lower boundary y=h(x0)
        
        if y > hx:

            return self.velocity(y, h0, hx)
        else: 
            return 0
        
    def velOutlet(self,y):
        h0 = self.yf           #upper boundary moving y=h0
        hx = self.hx_out        #lower boundary y=h(x0)
        
        if y > hx:
            return self.velocity(y, h0, hx)
        else: 

            return 0