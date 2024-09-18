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
    def __init__ (self, L, H, h, U, Q, Re, N, filestr):

        x0 = 0
        xf = L
        y0 = 0
        yf = H
        x_peaks = [x0, 1, xf]
        y_peaks=[[yf,yf-h],[yf-h,0],[0,yf]]

        super().__init__(x0, xf, y0, yf, N, U, Q, Re, filestr, x_peaks, y_peaks)
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
        Q = 2*U*h/2
        filestr = "stokes_BFS_H%dL%d_Re1_Q2_N%d"%(H, L, N)
        super().__init__(L, H, h, U, Q, Re, N, filestr)
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
        Q = 2*U*h/2
        filestr = "stokes_BFS_H%dL%d_Re05_N%d"%(H,L,N)
        super().__init__(L, H, h, U, Q, Re, N, filestr)
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
        Q = U*h/2 # --> dp/dx = 0 at inlet 
        filestr = "stokes_BFS_H%dL%d_Re0_Q1_N%d"%(H,L, N)
        super().__init__(L, H, h, U, Q, Re, N, filestr)

class BFS_H2L4_Re0_Q2(BFS):
    def __init__(self, N):
        L = 4
        H = 2
        h = 1   
        U = 1
        Re = 0
        Q = 2*U*h/2 
        filestr = "stokes_BFS_H%dL%d_Re0_Q2_N%d"%(H,L, N)
        super().__init__(L, H, h, U, Q, Re, N, filestr)
        
class BFS_H2L4_Re0_Q3(BFS):
    def __init__(self, N):
        L = 4
        H = 2
        h = 1   
        U = 1
        Re = 0
        Q = 3*U*h/2
        filestr = "stokes_BFS_H%dL%d_Re0_Q3_N%d"%(H,L, N)
        super().__init__(L, H, h, U, Q, Re, N, filestr)
        
#------------------------------------------------------------------------------
# <<====== Rect slider =======>>
#------------------------------------------------------------------------------

class RectSlider(PWLinear): 
    def __init__(self, L, l1, l2, H, h, U, Q, Re, N, filestr):
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

        super().__init__(x0, xf, y0, yf, N, U, Q, Re, filestr, x_peaks, y_peaks)


        
class RectSlider_H2L4_Re0_Q2(RectSlider):
    def __init__(self, N):
        L = 4
        l1 = 2
        l2 = 2
        H = 2
        h = 1
        U = 1
        Q = U*h
        Re = 0
        filestr = "stokes_rectSlider_rect_Re0_N%d"%N
        super().__init__(L, l1, l2, H, h, U, Q, Re, N, filestr)   





# class TriCavity(PWLinear):
#     def  __init__(self,L, H, U, Q, Re, N):
#         x0 = 0
#         xf = L
#         y0 = 0
#         yf = H
#         x_peaks = [x0, x0 + L/2, xf]
#         y_peaks = [[yf,yf],[y0,y0],[yf,yf]]
#         filestr = "stokes_slider_tri_Re1_N%d"%N
#         super().__init__(x0, xf, y0, yf, N, U, Q, Re, filestr, x_peaks, y_peaks)
        
        
        
# class Tri_standard(TriCavity):
#     def __init__(self,N):
#         Re = 1
#         U = 1
#         Q = 0
#         L = 1
#         H = 2
#         filestr = "stokes_slider_rectTxt_Re1_N%d"%N
#         super().__init__(L, l1, l2, H, h, U, Q, Re, N, filestr)

        
class trapSlider_Re0_Q2(PWLinear):
    def __init__(self, N):
        x0 = 0
        xf = 4
        y0 = 0
        yf=2
        x_peaks=[x0,1,2,3,xf]
        y_peaks = [[yf,1],[1,0.5],[0,0],[0.5,1],[1,yf]]
        h = 1
        U = 1
        Q = U*h
        Re = 0
        filestr = "stokes_slider_trap_Re0_N%d"%N
        super().__init__(x0, xf, y0, yf, N, U, Q, Re, filestr, x_peaks, y_peaks)

# class slider_triTxt_Re1(PWLinear):
#     def  __init__(self,N):
#         x0 = 0
#         xf = 10
#         y0 = 0
#         yf = 2
#         U = 1
#         Q = 0 
#         Re = 1 
#         x_peaks = [x0, 2, 3, 4, 5, 6, 7, xf]
#         y_peaks=[[yf,1],[1,1],[0,0],[1,1],[1,1],[0,0],[1,1],[1,yf]]
#         filestr = "stokes_slider_triTxt_Re1_N%d"%N
#         super().__init__(x0, xf, y0, yf, N, U, Q, Re, filestr, x_peaks, y_peaks)


# class slider_triTxt_Re05(PWLinear):
#     def  __init__(self,N):
#         x0 = 0
#         xf = 4
#         y0 = 0
#         yf = 2
#         U = 1
#         Q = 1
#         Re = 0.5 
#         x_peaks = [x0, 1, 2, 3, xf]
#         y_peaks=[[yf,1],[1,1],[0,0],[1,1],[1,yf]]
#         filestr = "stokes_slider_triTxt_Re05_N%d"%N
#         super().__init__(x0, xf, y0, yf, N, U, Q, Re, filestr, x_peaks, y_peaks)
  
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    