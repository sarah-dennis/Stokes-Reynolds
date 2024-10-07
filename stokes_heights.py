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
    
    def __init__(self, x0, xf, y0, yf, N, U, flux, Re, p0, filestr, x_peaks, y_peaks):
        super().__init__(x0, xf, y0, yf, N, U, flux, Re, p0, filestr)

        self.x_peaks = x_peaks
        self.y_peaks = y_peaks
        self.N_regions = len(self.x_peaks)-1
            
        # peaks must fall on the grid
        self.make_space()
        self.spacestr = "Textured Slider $Re=%.2f$"%Re  
                            
        # constants BCs on velocity, stream, flux 
        self.hf_in = y_peaks[0][1] # hf < h0 measured from y0
        self.hf_out = y_peaks[-1][0]
        
        self.H_in = yf - self.hf_in
        self.H_out = yf - self.hf_out
        
        if self.H_in == 0: #and self.H_out == 0: # closed cavity Q=0
            self.dp_in = 0
            # self.dp_out = 0 
        else:
            self.dp_in = (self.flux - 0.5*self.U*self.H_in) * (-12 / self.H_in**3)

            # self.dp_out = (self.flux - 0.5*self.U*self.H_out) * (-12 / self.H_out**3)
    
    # 0 : boundary, -1: exterior, 1: interior
    def make_space(self):

        slopes = np.zeros(self.N_regions)

        for k in range(self.N_regions):
            dh = self.y_peaks[k+1][0] - self.y_peaks[k][1]
            dx = self.x_peaks[k+1] - self.x_peaks[k]
            slopes[k] = dh/dx
    
        hs = np.zeros((self.Nx,2))
        
        grid = np.zeros((self.Ny, self.Nx))
    
        reg = 0 #pwl region        
        i_ref = 0
        for i in range(self.Nx):
            if self.xs[i] == self.x_peaks[reg]:
                h_left = self.y_peaks[reg][0]
                h_right = self.y_peaks[reg][1] 
                i_ref = i
                reg +=1
                hs[i] = [h_left,h_right] # choice for 1D
            else:
                h = slopes[reg-1]*(i - i_ref)*self.dx + self.y_peaks[reg-1][1]
                hs[i] = [h,h]
            
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
    
    # def streamOutlet(self, j):
    #     y = self.y0 + j*self.dy
    #     if y >= self.y_peaks[-1][0]:
    #         y = self.y0 + j*self.dy
    #         u_term = self.U* (0.5*(y**2 - self.yf**2) - self.hf_out*(y-self.yf))/self.H_out
    #         dp_term = -0.5*self.dp_out*( (-1/3)*(y**3 -self.yf**3) + 0.5*(self.yf+self.hf_out)*(y**2-self.yf**2) - self.yf*self.hf_out*(y-self.yf))
    #         return u_term + dp_term + self.flux
    #     else:
    #         return 0
        
    
    def velInlet(self,j):
        y = self.y0 + j*self.dy
        if y >= self.y_peaks[0][1]:
            
            u = (self.U/self.H_in - 0.5*self.dp_in*(self.yf-y)) * (y-self.hf_in) 

            return  u
        else: 
            return 0
        
    # def velOutlet(self,j):
    #     y = self.y0 + j*self.dy
    #     if y >= self.y_peaks[-1][0]:

    #         u = (self.U/self.H_out - 0.5*self.dp_out*(self.yf-y)) * (y-self.hf_out) 
    #         return u
    #     else: 
    #         return 0
#------------------------------------------------------------------------------
# Boundary interpolation
#------------------------------------------------------------------------------
 
# x_ij interior node
# x_st exterior nbr node 
    def interp_S(self, i,j, s,t, v_ij, v_bdry=0):
        y_ij = self.ys[j]
        y_nbr = self.ys[t]
        y_bdry = self.hs[i][0]
        
        l1 = (y_nbr-y_bdry)
        l2= -(y_ij-y_bdry) 
        v_nbr = v_bdry + (v_ij-v_bdry) * l1/l2

        return v_nbr

    def interp_E_W(self, i,j, s,t, v_ij, v_bdry=0):

        x_ij = self.xs[i]
        h_ij = self.hs[i][0]
        
        x_nbr = self.xs[s]
        if s == i + 1: #east
            h_nbr = self.hs[s][0]
        else: #west
            h_nbr = self.hs[s][1]

        y_bdry = self.ys[j]
        x_bdry = x_ij + (x_nbr-x_ij) * (y_bdry-h_ij)/(h_nbr-h_ij)
        l1 = (x_nbr-x_bdry)
        l2= -(x_ij-x_bdry)
        v_nbr = v_bdry +(v_ij-v_bdry) * l1/l2

        return v_nbr

    def interp_NE_SW(self, i,j, s,t, v_ij, v_bdry=0): 
        x_ij = self.xs[i]
        y_ij = self.ys[j]
        h_ij = self.hs[i][0]
        
        x_nbr = self.xs[s]
        y_nbr = self.ys[t]
        
        # if x_nbr, y_nbr  under a x=h(y) edge... 
        if s == i + 1: #east
            h_nbr = self.hs[s][0]
        else: #west
            h_nbr = self.hs[s][1]     
        slope = (h_nbr-h_ij)/(x_nbr-x_ij)
        x_bdry = (y_ij-h_ij)/(slope-1) + x_ij
        y_bdry = (x_bdry-x_ij) + y_ij
        
        l1 = np.sqrt((x_nbr-x_bdry)**2 + (y_nbr-y_bdry)**2)
        l2 = -np.sqrt((x_ij-x_bdry)**2 + (y_ij-y_bdry)**2)
        v_nbr = v_bdry + (v_ij-v_bdry)*l1/l2
        # print('nesw',x_bdry,y_bdry,v_nbr,v_ij)

        return v_nbr

    
    def interp_NW_SE(self, i,j, s,t, v_ij, v_bdry=0): 
        x_ij = self.xs[i]
        y_ij = self.ys[j]
        h_ij = self.hs[i][0]
        
        x_nbr = self.xs[s]
        y_nbr = self.ys[t]
        if s == i + 1: #east
            h_nbr = self.hs[s][0]
        else: #west
            h_nbr = self.hs[s][1]
        
        slope = (h_nbr-h_ij)/(x_nbr-x_ij)
        x_bdry = (y_ij-h_ij)/(slope+1) + x_ij
        y_bdry = -(x_bdry-x_ij) + y_ij
        
        l1 = np.sqrt((x_nbr-x_bdry)**2 + (y_nbr-y_bdry)**2)
        l2 = -np.sqrt((x_ij-x_bdry)**2 + (y_ij-y_bdry)**2)


        v_nbr = v_bdry + (v_ij-v_bdry)*l1/l2
        # print('nwse',x_bdry,y_bdry,v_nbr,v_ij)
        return v_nbr

#------------------------------------------------------------------------------
    def get_flux(self, delta_p):
        int_h_sqr,int_h_cube = self.integrate_h_flux()
        flux = ((-1/12)*delta_p - (1/2)*self.U*int_h_sqr )/int_h_cube
        return flux
    
    def integrate_h_flux(self):
        int_h_sqr = 0
        int_h_cube = 0
        for i in range(self.Nx-1):
            
            hi = self.hs[i][1]
            hj = self.hs[i+1][0]
            dh = hj - hi
            
            a = self.yf - 2*hi + hj
            b = self.yf - hi
            if dh == 0:
                int_h_cube += self.dx/b**3
                int_h_sqr += self.dx/b**2
            else:
                int_h_cube += -0.5*self.dx/dh * (1/a**2 - 1/b**2) 
                int_h_sqr += -self.dx/dh *(1/a - 1/b) 
        return int_h_sqr, int_h_cube


      























