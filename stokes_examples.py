# -*- coding: utf-8 -*-
"""
Created on Tue May 21 16:15:37 2024

@author: sarah
"""

from stokes_heights import triangle
from stokes_heights import step
# example = examples.tri_Re1
# example = examples.tri_Re0

# example = examples.bfs_Re1
# example = examples.bfs_Re10neg4


class tri_Re1(triangle):
    def __init__(self, N):
        x0 = 0
        xf = 1
        y0 = 0
        yf = 2
        U = 1 # velocity @ yf
        Q = 0 # stream @ yf
        Re = 1
        filestr = "stokes_tri_Re1_N%d"%(N)
        super().__init__(x0, xf, y0, yf, N, U, Q, Re, filestr)
        
class tri_Re0(triangle):
    def __init__(self, N):
        x0 = 0
        xf = 1
        y0 = 0
        yf = 2
        U = 1
        Q = 0
        Re = 0
        filestr = "stokes_tri_Re0_N%d"%(N)
        super().__init__(x0, xf, y0, yf, N, U, Q, Re, filestr)
        
class bfs_Re10neg4(step):
    def __init__(self, N):
        x0 = 0
        xf = 4
        y0 = 0
        yf = 2
        U = 0
        Q = 1e-4 # stream @ yf
        Re = 1e-4
        x_step = 4
        y_step = 2
        filestr = "stokes_BFS_Re1e-4_N%d"%N

        super().__init__(x0, xf, y0, yf, N, U, Q, Re, filestr, x_step, y_step)
        
class bfs_ReHalf(step):
    def __init__(self, N):
        x0 = 0
        xf = 4
        y0 = 0
        yf = 2
        U = 1
        Re = 0.5
        Q = 0 # * (viscosity / density)
        x_step = 4
        y_step = 2
        filestr = "stokes_BFS_Refe-1_N%d"%N

        super().__init__(x0, xf, y0, yf, N, U, Q, Re, filestr, x_step, y_step)