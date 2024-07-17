# -*- coding: utf-8 -*-
"""
Created on Tue May 21 16:15:37 2024

@author: sarah
"""

from stokes_heights import triangle
from stokes_heights import step

class biswasEx(triangle):
    def __init__(self, N):
        x0 = 0
        xf = 1
        y0 = 0
        yf = 2
        U = 1
        Re = 1
        filestr = "stokes_tri_N%d_spLU"%(N)
        super().__init__(x0, xf, y0, yf, N, U, Re, filestr)
        
class zeroReynEx(triangle):
    def __init__(self, N):
        x0 = 0
        xf = 1
        y0 = 0
        yf = 2
        U = 1
        Re = 0
        filestr = "stokes_tri2_N%d_spLU"%(N)
        super().__init__(x0, xf, y0, yf, N, U, Re, filestr)
        
        
class bfsEx1(step):
    def __init__(self, N):
        x0 = 0
        xf = 10
        y0 = 0
        yf = 1
        U = 0
        Q = 1
        Re = 0.1
        x_step = 2
        y_step = 2
        filestr = "stokes_BFS1_N%d_spLU"%(N)
        super().__init__(x0, xf, y0, yf, N, U, Q, Re, filestr, x_step, y_step)
        
class bfsEx2(step):
    def __init__(self, N):
        x0 = 0
        xf = 10
        y0 = 0
        yf = 1
        U = 0
        Q = 10
        Re = 0.1
        x_step = 2
        y_step = 2
        filestr = "stokes_BFS2_N%d_spLU"%(N)
        super().__init__(x0, xf, y0, yf, N, U, Q, Re, filestr, x_step, y_step)

        
class bfsEx3(step):
    def __init__(self, N):
        x0 = 0
        xf = 10
        y0 = 0
        yf = 1
        U = 5
        Q = 10
        Re = 0.1
        x_step = 5
        y_step = 2
        filestr = "stokes_BFS3_N%d_spLU"%(N)
        super().__init__(x0, xf, y0, yf, N, U, Q, Re, filestr, x_step, y_step)
    
class bfs_biswasLowRe(step):
    def __init__(self, N):
        x0 = 0
        xf = 4
        y0 = 0
        yf = 2
        U = 0
        Re = 1e-2
        Q = Re # * (viscosity / density)
        x_step = 4
        y_step = 2
        filestr = "stokes_BFS4_N%d_spLU"%N

        super().__init__(x0, xf, y0, yf, N, U, Q, Re, filestr, x_step, y_step)
        
class bfs_biswasLowerRe(step):
    def __init__(self, N):
        x0 = 0
        xf = 4
        y0 = 0
        yf = 2
        U = 0
        Re = 1e-4
        Q = Re # * (viscosity / density)
        x_step = 4
        y_step = 2
        filestr = "stokes_BFS5_N%d_spLU"%N

        super().__init__(x0, xf, y0, yf, N, U, Q, Re, filestr, x_step, y_step)
        
        
