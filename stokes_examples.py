 # -*- coding: utf-8 -*-
"""
Created on Tue May 21 16:15:37 2024

@author: sarah
"""

from stokes_heights import PWLinear
import numpy as np
#------------------------------------------------------------------------------
# <<===== BFS ======>>
#------------------------------------------------------------------------------
class BFS(PWLinear):
    def __init__ (self, L, H, h, U, q, Re, N, filestr):
        p_amb = 0
        x0 = 0
        xf = L
        y0 = 0
        yf = H
        x_peaks = [x0, 1, xf]
        y_peaks=[[yf,yf-h],[yf-h,0],[0,yf]]
        super().__init__(x0, xf, y0, yf, N, U, q, Re, p_amb, filestr, x_peaks, y_peaks)
#------------------------------------------------------------------------------
# BFS Re=1
#------------------------------------------------------------------------------

class BFS_H2L4_Re1_Q2(BFS):
    def __init__(self, N):
        L = 4
        H = 2
        h = 1   
        U = 1
        Re = 1
        q = 2
        filestr = "stokes_BFS_H%dL%d_Re1_Q2_N%d"%(H, L, N)
        super().__init__(L, H, h, U, q, Re, N, filestr)
#------------------------------------------------------------------------------
# BFS Re=1/2
#------------------------------------------------------------------------------
class BFS_H2L4_Re05_Q2(BFS):
    def __init__(self, N):
        L = 4
        H = 2
        h = 1   
        U = 1
        Re = 0.5
        q = 2
        filestr = "./examples/stokes_BFS_H%dL%d_Re05_Q2_N%d"%(H,L,N)
        super().__init__(L, H, h, U, q, Re, N, filestr)
#------------------------------------------------------------------------------
# BFS Re=0
#------------------------------------------------------------------------------       
class BFS_H2L4_Re0_Q1(BFS):
    def __init__(self, N):
        L = 4
        H = 2
        h = 1   
        U = 1
        Re = 0
        q = 1 
        filestr = "./examples/stokes_BFS_H%dL%d_Re0_Q1_N%d"%(H,L, N)
        super().__init__(L, H, h, U, q, Re, N, filestr)

class BFS_H2L4_Re0_Q2(BFS):
    def __init__(self, N):
        L = 4
        H = 2
        h = 1   
        U = 1
        Re = 0
        q = 2
        filestr = "stokes_BFS_H%dL%d_Re0_Q2_N%d"%(H,L, N)
        super().__init__(L, H, h, U, q, Re, N, filestr)
        

        
#------------------------------------------------------------------------------
# <<====== Rect slider =======>>
#------------------------------------------------------------------------------

class RectSlider(PWLinear): 
    def __init__(self, L, l1, l2, H, h, U, q, Re, N, filestr):
        p_amb=0
        x0 = 0
        xf = L #+x0
        y0 = 0
        yf = H #+y0
        N_texts = L//(l1 + l2)
        x_peaks = np.zeros(2*N_texts+2)
        y_peaks = np.zeros((2*N_texts+2, 2))      
        x_peaks[0] = x0
        y_peaks[0] = [yf,yf-h]
        x_peaks[1] = l1/2
        y_peaks[1] = [yf-h, y0]
        x_peaks[-1] = xf
        y_peaks[-1] = [yf-h, yf]
        for i in range(2, 2*N_texts+1):
            if i % 2 == 1:
                x_peaks[i] = x_peaks[i-1] + l1
                y_peaks[i] = [yf-h, y0]
            else:
                x_peaks[i] = x_peaks[i-1] + l2
                y_peaks[i] = [y0, yf-h]

        super().__init__(x0, xf, y0, yf, N, U, q, Re, p_amb, filestr, x_peaks, y_peaks)


        
class RectSlider_H2L4_Re0_Q2(RectSlider):
    def __init__(self, N):
        L = 4
        l1 = 2
        l2 = 2
        H = 2
        h = 1
        U = 1
        q = 2
        Re = 0
        filestr = "stokes_slider_rect_Re0_N%d"%N
        super().__init__(L, l1, l2, H, h, U, q, Re, N, filestr)   
        
class RectSlider_H2L4_Re05_Q2(RectSlider):
    def __init__(self, N):
        L = 4
        l1 = 2
        l2 = 2
        H = 2
        h = 1
        U = 1
        q = 2
        Re = 0.5
        filestr = "stokes_slider_rect_Re0_N%d"%N
        super().__init__(L, l1, l2, H, h, U, q, Re, N, filestr)   

#------------------------------------------------------------------------------
# <<====== Trapezoid slider =======>>
#------------------------------------------------------------------------------

class HexSlider_Re05_Q2(PWLinear):
    def __init__(self, N):
        x0 = 0
        xf = 4
        y0 = 0
        yf = 2
        x_peaks = [x0,1,2,3,xf]
        y_peaks = [[yf,1],[1,0.5],[0,0],[0.5,1],[1,yf]]
        U = 0
        q = 2
        p0 = 10
        Re = 0.5
        filestr = "./examples/stokes_slider_hex_Re05_Q2_N%d"%N
        super().__init__(x0, xf, y0, yf, N, U, q, Re, p0, filestr, x_peaks, y_peaks)
      

#------------------------------------------------------------------------------
# <<====== Triangle slider =======>>
#------------------------------------------------------------------------------
class TriSlider_Re05_Q1(PWLinear):
    def  __init__(self,N):
        L=4
        l=1
        H=2
        h=1
        x0 = 0
        xf = L
        y0 = 0
        yf = H
        Re=0.5
        U=0.5 
        q=1
        p0=0
        x_peaks = [x0, x0 + l, x0 + 2*l, x0 +3*l, xf]
        y_peaks = [[yf,yf-h],[yf-h,yf-h],[y0,y0],[yf-h,yf-h],[yf-h,yf]]
        filestr = "stokes_slider_tri_Re0_Q2_N%d"%N
        super().__init__(x0, xf, y0, yf, N, U, q, Re, p0, filestr, x_peaks, y_peaks)
        
        
    
    
    
    
    
    
    
    
    
    
    