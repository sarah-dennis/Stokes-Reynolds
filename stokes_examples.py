# -*- coding: utf-8 -*-
"""
Created on Tue May 21 16:15:37 2024

@author: sarah
"""

from stokes_heights import triangle

class biswasEx(triangle):
    def __init__(self, N):
        x0 = 0
        xf = 1
        y0 = 0
        yf = 2
        U = 1
        Re = 1
        filestr = "stokes_N%d_spLU"%(N)
        super().__init__(x0, xf, y0, yf, N, U, Re, filestr)
        
class zeroReynEx(triangle):
    def __init__(self, N):
        x0 = 0
        xf = 1
        y0 = 0
        yf = 2
        U = 1
        Re = 0
        filestr = "stokes_Re0_N%d_spLU"%(N)
        super().__init__(x0, xf, y0, yf, N, U, Re, filestr)
        