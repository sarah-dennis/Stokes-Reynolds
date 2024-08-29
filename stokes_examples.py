 # -*- coding: utf-8 -*-
"""
Created on Tue May 21 16:15:37 2024

@author: sarah
"""

from stokes_heights import PWLinear
import numpy as np
#------------------------------------------------------------------------------
class BFS(PWLinear):
    def __init__ (self, L, H, h, U, Q, Re, N):

        x0 = 0
        xf = L
        y0 = 0
        yf = H
        x_peaks = [x0, L/5, xf]
        y_peaks=[[yf,yf-h],[yf-h,0],[0,yf]]
        filestr = "stokes_BFS_H%d_Re%d_N%d"%(H, Re, N)
        super().__init__(x0, xf, y0, yf, N, U, Q, Re, filestr, x_peaks, y_peaks)

class BFS_standard(BFS):
    def __init__(self, N):
        Re = 0
        U = 1
        Q = 1
        L = 5
        H = 2
        h = 1
        super().__init__(L, H, h, U, Q, Re, N)

#------------------------------------------------------------------------------

class TriCavity(PWLinear):
    def  __init__(self,L, H, U, Q, Re, N):
        x0 = 0
        xf = L
        y0 = 0
        yf = H
        x_peaks = [x0, L/2, xf]#, 3, 4,5]
        y_peaks = [[yf,H],[0,0],[H,yf]]
        filestr = "stokes_slider_tri_Re1_N%d"%N
        super().__init__(x0, xf, y0, yf, N, U, Q, Re, filestr, x_peaks, y_peaks)
        
class Tri_standard(TriCavity):
    def __init__(self,N):
        Re = 1
        U = 1
        Q = 0
        L = 1
        H = 2
        super().__init__(L, H, U, Q, Re, N)

#------------------------------------------------------------------------------
class RectTexturedSlider(PWLinear):
    def __init__(self, L, l1, l2, H, h, U, Q, Re, N):
        x0 = 0
        xf = L
        y0 = 0
        yf = H
        N_texts = L//(l1 + l2)
        x_peaks = np.zeros(2*N_texts+2)
        y_peaks = np.zeros((2*N_texts+2, 2))      
        x_peaks[0] = x0
        y_peaks[0] = [yf,H-h]
        x_peaks[1] = l1/2
        y_peaks[1] = [H-h, y0]
        x_peaks[-1] = xf
        y_peaks[-1] = [H-h, yf]
        for i in range(2, 2*N_texts+1):
            if i % 2 == 1:
                x_peaks[i] = x_peaks[i-1] + l1
                y_peaks[i] = [H-h, y0]
            else:
                x_peaks[i] = x_peaks[i-1] + l2
                y_peaks[i] = [y0, H-h]
        filestr = "stokes_slider_rectTxt_Re1_N%d"%N
        super().__init__(x0, xf, y0, yf, N, U, Q, Re, filestr, x_peaks, y_peaks)


class rectSlider_standard(RectTexturedSlider):
    def __init__(self, N):
        L = 20
        l1 = 2
        l2 = 3
        H = 2
        h = 1
        U = 1
        Q = 1
        Re = 1
        super().__init__(L, l1, l2, H, h, U, Q, Re, N)
    
class slider_rectTxt_Re1(PWLinear):
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

class slider_triTxt_Re1(PWLinear):
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


    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    