 # -*- coding: utf-8 -*-
"""
Created on Tue May 21 16:15:37 2024

@author: sarah
"""

from stokes_heights import triangle
from stokes_heights import step
from stokes_heights import slider

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
        
        
#------------------------------------------------------------------------------

class slider_Re1(slider):
    def  __init__(self,N):
        x0 = 0
        xf = 2
        y0 = 0
        yf = 2
        U = 0
        Q = 1 
        Re = 1 
        x_peaks = [x0, 1, xf]#, 3, 4,5]
        y_peaks=[[yf,1],[0.5,1],[3/2,yf]]#,[0.5,0.5],[0.25,0.5],[1,0]]
        filestr = "stokes_slider_Re1_N%d"%N
        super().__init__(x0, xf, y0, yf, N, U, Q, Re, filestr, x_peaks, y_peaks)

class slider_tri_Re1(slider):
    def  __init__(self,N):
        x0 = 0
        xf = 1
        y0 = 0
        yf = 2
        U = 1
        Q = 0 
        Re = 1 
        x_peaks = [x0, 0.5, xf]#, 3, 4,5]
        y_peaks=[[yf,2],[0,0],[2,yf]]
        filestr = "stokes_slider_tri_Re1_N%d"%N
        super().__init__(x0, xf, y0, yf, N, U, Q, Re, filestr, x_peaks, y_peaks)

class slider_rectTxt_Re1(slider):
    def  __init__(self,N):
        x0 = 0
        xf = 10
        y0 = 0
        yf = 2
        U = 1
        Q = 0 
        Re = 1 
        x_peaks = [x0, 2, 2.5, 3.5, 4, 5, 5.5, 6.5, 7, xf]#, 3, 4,5]
        y_peaks=[[yf,1],[1,0],[0,1],[1,0],[0,1],[1,0],[0,1],[1,0],[0,1],[1,yf]]
        filestr = "stokes_slider_rectTxt_Re1_N%d"%N
        super().__init__(x0, xf, y0, yf, N, U, Q, Re, filestr, x_peaks, y_peaks)

class slider_triTxt_Re1(slider):
    def  __init__(self,N):
        x0 = 0
        xf = 10
        y0 = 0
        yf = 2
        U = 1
        Q = 0 
        Re = 1 
        x_peaks = [x0, 2, 3, 4, 5, 6, 7, xf]#, 3, 4,5]
        y_peaks=[[yf,1],[1,1],[0,0],[1,1],[1,1],[0,0],[1,1],[1,yf]]
        filestr = "stokes_slider_triTxt_Re1_N%d"%N
        super().__init__(x0, xf, y0, yf, N, U, Q, Re, filestr, x_peaks, y_peaks)


    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    