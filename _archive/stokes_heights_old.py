# -*- coding: utf-8 -*-
"""
Created on Tue Aug 27 08:25:06 2024

@author: sarah
"""

class tri_Re1(triangle):
    def __init__(self, N):
        x0 = 0
        xf = 1
        y0 = 0
        yf = 2
        slopes = [-4,4]
        wavelen = 1 * (xf-x0)
        U = 1 # velocity @ yf
        Q = 0 # stream @ yf
        Re = 1 # rhs factor
        filestr = "stokes_tri_Re1_N%d"%(N)
        super().__init__(x0, xf, y0, yf, N, U, Q, Re, filestr, slopes, wavelen)
        
class tri_Re0(triangle):
    def __init__(self, N):
        x0 = 0
        xf = 1                       
        y0 = 0
        yf = 2
        slopes = [-4,4]
        wavelen = 1 * (xf-x0)
        U = 1 # velocity @ yf
        Q = 0
        Re = 0
        filestr = "stokes_tri_Re0_N%d"%(N)
        super().__init__(x0, xf, y0, yf, N, U, Q, Re, filestr, slopes, wavelen)
#------------------------------------------------------------------------------
class bfs_Re1(step):
    def __init__(self, N):
        x0 = 0
        xf = 4
        y0 = 0
        yf = 2
        U = 0 # velocity @ yf
        Q = 1 # stream @ yf
        Re = 1
        x_step = 4 #Lx_out/Lx_in
        y_step = 2 #Ly_out/Ly_in
        filestr = "stokes_BFS_Re1_N%d"%N
        super().__init__(x0, xf, y0, yf, N, U, Q, Re, filestr, x_step, y_step)
        
class bfs_Re0(step):
    def __init__(self, N):
        x0 = 0
        xf = 4
        y0 = 0
        yf = 2
        U = 0 # velocity @ yf
        Q = 1 # stream @ yf
        Re = 0
        x_step = 4
        y_step = 2
        filestr = "stokes_BFS_Re1_N%d"%N
        super().__init__(x0, xf, y0, yf, N, U, Q, Re, filestr, x_step, y_step)
                

class bfs_Re10neg4(step):
    def __init__(self, N):
        x0 = 0
        xf = 4
        y0 = 0
        yf = 2
        U = 0 # velocity @ yf
        Q = 1e-4 # stream @ yf
        Re = 1e-4
        x_step = 4
        y_step = 2
        filestr = "stokes_BFS_Re1e-4_N%d"%N
        super().__init__(x0, xf, y0, yf, N, U, Q, Re, filestr, x_step, y_step)
        
        
class triangle(Space):
    def __init__(self, x0, xf, y0, yf, N, U, Q, Re, filestr, slopes, wavelen):

        super().__init__(x0, xf, y0, yf, N, U, Q, Re, filestr)
        self.slope_A = slopes[0]
        self.slope_B = slopes[1]
        self.N_regions = 2
        self.slopes = slopes
        self.x_peaks = [x0, (xf-x0)/(2*wavelen), xf]
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
        # graphics.plot_contour_mesh(grid, self.xs, self.ys, 'space',['space', 'x', 'y'])

        return grid
    
    
    def interp_E_W(self, t, s, psi_k):
        x = self.xs[s]
        for k in range(self.N_regions):
            if x >= self.x_peaks[k]:
                slope = self.slopes[k]
                break
        
        scale = 1 - (t%slope)/slope
        
        return -scale * psi_k

    def interp_S(self, t, s, psi_k):
        x = self.xs[s]
        for k in range(self.N_regions):
            if x >= self.x_peaks[k]:
                slope = self.slopes[k]
                break
        scale = 1 - (s % (1/slope))*slope
        return -scale * psi_k

    def interp_NE_NW(self, t, s, psi_N):
        x = self.xs[s]
        for k in range(self.N_regions):
            if x >= self.x_peaks[k]:
                slope = self.slopes[k]
                break
        
        scale = 1 - (t%slope)/slope
        
        return -scale * psi_N

    def interp_SE_SW(self, t, s, psi_SEW, psi_S):
        x = self.xs[s]
        for k in range(self.N_regions):
            if x >= self.x_peaks[k]:
                slope = self.slopes[k]
                break
        
        if abs(slope) < 1:
            
            scale = 1 - (s%(1/slope))*slope
        
            return -scale * psi_SEW
        else:
            scale = 1 - (t % slope)/slope
            return -scale * psi_S

    

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
        