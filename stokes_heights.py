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
    
    def __init__(self, x0, xf, y0, yf, N, U, Q, Re, p_amb, filestr, x_peaks, y_peaks):
        super().__init__(x0, xf, y0, yf, N, U, Q, Re, p_amb, filestr)
        # peaks must fall on the grid
        self.x_peaks = x_peaks
        self.y_peaks = y_peaks
        self.N_regions = len(self.x_peaks)-1
        self.make_space()
        self.spacestr = "$Re=%.2f$, $Q=%.2f$, $U=%.1f$"%(Re,Q,U)  
                            
        self.hf_in = y_peaks[0][1] # hf < h0 measured from y0
        self.H_in = yf - self.hf_in
        self.H_out = yf-y_peaks[-1][0]
        
        if self.H_in == 0: # closed cavity --> Q=0, dp=0
            self.dp_in = 0
        else: # gap entry --> dp ~ Q
            self.dp_in = (self.flux - 0.5*self.U*self.H_in) * (-12 / self.H_in**3)
    
    # 0 : boundary, -1: exterior, 1: interior
    def make_space(self):

        slopes = np.zeros(self.N_regions)

        for k in range(self.N_regions):
            dh = self.y_peaks[k+1][0] - self.y_peaks[k][1]
            dx = self.x_peaks[k+1] - self.x_peaks[k]
            slopes[k] = dh/dx
        hs = np.zeros((self.Nx,2))
        
        grid = np.zeros((self.Ny, self.Nx))
    
        reg = 0        
        i_ref = 0
        for i in range(self.Nx):
            if math.isclose(self.xs[i], self.x_peaks[reg]):
                h_left = self.y_peaks[reg][0]
                h_right = self.y_peaks[reg][1] 
                i_ref = i
                reg +=1
                hs[i] = [h_left,h_right]
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
                    if math.isclose(y, h_left) or math.isclose(y, h_right): # true boundary point at region change
                        grid[j,i] = 0
                        
                    elif h_left < h_right:
                        if h_left < y and y < h_right: # x=h(y) vertical boundary
                            grid[j,i] = 0
                        elif y > h_right: # above vert boundary (interior)
                            grid[j,i] = 1
                        else:             # below vert boundary(exterior)
                            grid[j,i] = -1
                
                    else:
                        if h_left > y and y > h_right: # x=h(y) vertical boundary
                            grid[j,i] = 0
                        elif y > h_left: # above vert boundary (interior)
                            grid[j,i] = 1
                        else:              # below vert boundary (exterior)
                            grid[j,i] = -1

                else:
                    if math.isclose(y,h): # true boundary point not at region change (from dx | slope)
                        grid[j,i] = 0
                    elif y > h:            # above boundary (interior)

                        grid[j,i] = 1
                    else:                   # below boundary (exterior)
                        grid[j,i] = -1
        
        self.space = grid
        self.slopes = slopes
        self.hs = hs

#------------------------------------------------------------------------------
# Boundary conditions on stream and velocity
#------------------------------------------------------------------------------
   
    def streamInlet(self, j):
        y = self.y0 + j*self.dy
        if y >= self.y_peaks[0][1]:
            u_term = self.U*(0.5*(y**2 - self.yf**2) - self.hf_in*(y-self.yf))/self.H_in
            dp_term = -0.5*self.dp_in*((-1/3)*(y**3 -self.yf**3) + 0.5*(self.yf+self.hf_in)*(y**2-self.yf**2) - self.yf*self.hf_in*(y-self.yf))
            psi = u_term + dp_term + self.flux 
            return psi
        else:
            return 0
    
    def velInlet(self, j):
        y = self.y0 + j*self.dy
        if y >= self.y_peaks[0][1]:
            u = (self.U/self.H_in - 0.5*self.dp_in*(self.yf-y)) * (y-self.hf_in) 
            return  u
        else: 
            return 0
        
#------------------------------------------------------------------------------
# Boundary interpolation
#------------------------------------------------------------------------------

# x_ij interior node
# x_st exterior nbr node 
    def interp_S(self, i,j, s,t, v_ij, v_bdry=0):
        if np.isclose(v_ij,v_bdry):
            return v_bdry
        
        y_ij = self.ys[j]
        y_nbr = self.ys[t]
        y_bdry = self.hs[i][0] #arbitrary -- South lookup not at x-peaks
        
        l1 = (y_nbr-y_bdry)
        l2 = (y_ij-y_bdry) 
        v_nbr = v_bdry + (v_ij-v_bdry) * (l1/l2)        
        assert v_nbr == 0 or np.sign(v_nbr) != np.sign(v_ij), "s"
        return v_nbr

    def interp_E_W(self, i,j, s,t, v_ij, v_bdry=0):
        if np.isclose(v_ij,v_bdry):
            return v_bdry
        
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
        l2 = (x_ij-x_bdry)
        v_nbr = v_bdry + (v_ij-v_bdry) * (l1/l2)
        assert v_nbr == 0 or np.sign(v_nbr) != np.sign(v_ij), "e-w"
        return v_nbr

    def interp_NE_SW(self, i,j, s,t, v_ij, v_bdry=0): 
        if np.isclose(v_ij,v_bdry):
            return v_bdry
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
        x_bdry = (y_ij-h_ij)/(slope-1) + x_ij
        y_bdry = (x_bdry-x_ij) + y_ij
        l1 = np.sqrt((x_nbr-x_bdry)**2 + (y_nbr-y_bdry)**2)
        l2 = -np.sqrt((x_ij-x_bdry)**2 + (y_ij-y_bdry)**2)
        
        v_nbr = v_bdry + (v_ij-v_bdry)* (l1/l2)
        assert v_nbr == 0 or np.sign(v_nbr) != np.sign(v_ij), "ne-sw"
        return v_nbr

    
    def interp_NW_SE(self, i,j, s,t, v_ij, v_bdry=0): 
        if np.isclose(v_ij,v_bdry):
            return v_bdry
        
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

        v_nbr = v_bdry + (v_ij-v_bdry)*(l1/l2)
        assert v_nbr == 0 or np.sign(v_nbr) != np.sign(v_ij), "nw-se"
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





















