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
    def __init__ (self, L, l, H, h, U, q, Re, N, filestr):
        p_amb = 0
        x0 = 0
        xf = L
        y0 = 0
        yf = H
        x_peaks = [x0, x0+l, xf]
        y_peaks=[[yf,yf-h],[yf-h,0],[0,yf]]
        super().__init__(x0, xf, y0, yf, N, U, q, Re, p_amb, filestr, x_peaks, y_peaks)
#------------------------------------------------------------------------------
        
class BFS_H2L4_Re1_Q2_U0(BFS):
    def __init__(self, N):
        L = 4
        l=1
        H = 2
        h = 1   
        U = 0
        Re = 1
        q = 2
        name_str = "BFS_H2L4_Re1_Q2_U0"
        filestr = "./examples/%s/%s_N%s"%(name_str, name_str, N)
        super().__init__(L, l, H, h, U, q, Re, N, filestr)
        
class BFS_H2L4_Re0p5_Q2_U0(BFS):
    def __init__(self, N):
        L = 4
        l=1
        H = 2
        h = 1   
        U = 0
        Re = 0.5
        q = 2
        name_str = "BFS_H2L4_Re0p5_Q2_U0"
        filestr = "./examples/%s/%s_N%s"%(name_str, name_str, N)
        super().__init__(L, l, H, h, U, q, Re, N, filestr)

        
class BFS_H2L4_Re0_Q2_U0(BFS):
    def __init__(self, N):
        L = 4
        l=1
        H = 2
        h = 1   
        U = 0
        Re = 0
        q = 2
        name_str = "BFS_H2L4_Re0_Q2_U0"
        filestr = "./examples/%s/%s_N%s"%(name_str, name_str, N)
        super().__init__(L, l, H, h, U, q, Re, N, filestr)
        
class BFS_H1p5L4_Re0_Q2_U0(BFS):
    def __init__(self, N):
        L = 4
        l=1
        H = 1.5
        h = 1   
        U = 0
        Re = 0
        q = 2
        name_str = "BFS_H1p5L4_Re0_Q2_U0"
        filestr = "./examples/%s/%s_N%s"%(name_str, name_str, N)
        super().__init__(L, l, H, h, U, q, Re, N, filestr)
        
class BFS_H1p25L4_Re0_Q2_U0(BFS):
    def __init__(self, N):
        L = 4
        l = 1
        H = 1.25
        h = 1   
        U = 0
        Re = 0
        q = 2
        name_str = "BFS_H1p25L4_Re0_Q2_U0"
        filestr = "./examples/%s/%s_N%s"%(name_str, name_str, N)
        super().__init__(L, l, H, h, U, q, Re, N, filestr)
         
class BFS_H1p125L4_Re0_Q2_U0(BFS):
    def __init__(self, N):
        L = 4
        l = 1
        H = 1.125
        h = 1   
        U = 0
        Re = 0
        q = 2
        name_str = "BFS_H1p125L4_Re0_Q2_U0"
        filestr = "./examples/%s/%s_N%s"%(name_str, name_str, N)
        super().__init__(L, l, H, h, U, q, Re, N, filestr)
#------------------------------------------------------------------------------
# <<====== Smoothed BFS =======>>
#------------------------------------------------------------------------------
class BFS_deltaSmooth(PWLinear):
    def __init__ (self, L, l,H, h, delta, U, q, Re, N, filestr):
        p_amb = 0
        x0 = 0
        xf = L
        y0 = 0
        yf = H

        x_peaks = [x0, l-delta, l, l+delta, xf]
        y_peaks=[[yf,yf-h],[yf-h,yf-h],[yf-h-(H-h)/2,yf-h-(H-h)/2],[0,0],[0,yf]]
        super().__init__(x0, xf, y0, yf, N, U, q, Re, p_amb, filestr, x_peaks, y_peaks)



class BFS_H2L4_delta1p5(BFS_deltaSmooth):
    def __init__ (self, N):
        L=4
        l=2
        H= 2 
        h=1 
        delta=1.5
        U=0
        q=2
        Re=0
        name_str = "dBFS_H2L4_d1p5_Re0_Q2_U0"
        filestr = "./examples/%s/%s_N%s"%(name_str, name_str, N)
        super().__init__(L, l, H, h, delta, U, q, Re, N, filestr)


class BFS_H2L4_delta1(BFS_deltaSmooth):
    def __init__ (self, N):
        L=4
        l=2
        H= 2 
        h=1 
        delta=1
        U=0
        q=2
        Re=0
        name_str = "dBFS_H2L4_d1_Re0_Q2_U0"
        filestr = "./examples/%s/%s_N%s"%(name_str, name_str, N)
        super().__init__(L, l, H, h, delta, U, q, Re, N, filestr)

class BFS_H2L4_delta0p75(BFS_deltaSmooth):
    def __init__ (self, N):
        L=4
        l=2
        H= 2 
        h=1 
        delta=0.75
        U=0
        q=2
        Re=0
        name_str = "dBFS_H2L4_d0p75_Re0_Q2_U0"
        filestr = "./examples/%s/%s_N%s"%(name_str, name_str, N)
        super().__init__(L, l, H, h, delta, U, q, Re, N, filestr)
        
class BFS_H2L4_delta0p5(BFS_deltaSmooth):
    def __init__ (self, N):
        L=4
        l=2
        H= 2 
        h=1 
        delta=0.5
        U=0
        q=2 
        Re=0
        name_str = "dBFS_H2L4_d0p5_Re0_Q2_U0"
        filestr = "./examples/%s/%s_N%s"%(name_str, name_str, N)
        super().__init__(L,l, H, h, delta, U, q, Re, N, filestr)

class BFS_H2L4_delta0p25(BFS_deltaSmooth):
    def __init__ (self, N):
        L=4
        l=2
        H= 2 
        h=1 
        delta=0.25
        U=0
        q=2 
        Re=0
        folder_str = "BFS_H2L4_Re0_Q2_U0_deltas/delta0p25"
        name_str = "BFS_H2L4_delta0p25_Re0_Q2_U0"
        filestr = "./examples/%s/%s_N%s"%(folder_str, name_str, N)
        super().__init__(L, l, H, h, delta, U, q, Re, N, filestr)
        
        #BLOWS UP!!
class BFS_H2L4_delta0p125(BFS_deltaSmooth):
    def __init__ (self, N):
        L=4
        l=1
        H= 2 
        h=1 
        delta=0.125
        U=0
        q=2 
        Re=0
        name_str = "dBFS_H2L4_d0p125_Re0_Q2_U0"
        filestr = "./examples/%s/%s_N%s"%(name_str, name_str, N)
        super().__init__(L, l, H, h, delta, U, q, Re, N, filestr)

class BFS_H2L4_delta0(BFS):
    def __init__ (self, N):
        L=4
        l=2
        H= 2 
        h=1 
        U=0
        q=2 
        Re=0
        name_str = "dBFS_H2L4_d0_Re0_Q2_U0"
        filestr = "./examples/%s/%s_N%s"%(name_str, name_str, N)
        super().__init__(L, l, H, h, U, q, Re, N, filestr)
        


#------------------------------------------------------------------------------



class BFS_H1p5L4_delta1p25(BFS_deltaSmooth):
    def __init__ (self, N):
        L=4
        l=2
        H=1.5
        h=1 
        delta=1.25
        U=0
        q=2
        Re=0
        name_str = "dBFS_H1p5L4_d1p25_Re0_Q2_U0"
        filestr = "./examples/%s/%s_N%s"%(name_str, name_str, N)
        super().__init__(L, l, H, h, delta, U, q, Re, N, filestr)


class BFS_H1p5L4_delta1(BFS_deltaSmooth):
    def __init__ (self, N):
        L=4
        l=2
        H= 1.5
        h=1 
        delta=1
        U=0
        q=2
        Re=0
        name_str = "dBFS_H1p5L4_d1_Re0_Q2_U0"
        filestr = "./examples/%s/%s_N%s"%(name_str, name_str, N)
        super().__init__(L, l, H, h, delta, U, q, Re, N, filestr)

class BFS_H1p5L4_delta0p75(BFS_deltaSmooth):
    def __init__ (self, N):
        L=4
        l=2
        H=1.5 
        h=1 
        delta=0.75
        U=0
        q=2
        Re=0
        name_str = "dBFS_H1p5L4_d0p75_Re0_Q2_U0"
        filestr = "./examples/%s/%s_N%s"%(name_str, name_str, N)
        super().__init__(L, l, H, h, delta, U, q, Re, N, filestr)
        
class BFS_H1p5L4_delta0p5(BFS_deltaSmooth):
    def __init__ (self, N):
        L=4
        l=2
        H=1.5
        h=1 
        delta=0.5
        U=0
        q=2 
        Re=0
        name_str = "dBFS_H1p5L4_d0p5_Re0_Q2_U0"
        filestr = "./examples/%s/%s_N%s"%(name_str, name_str, N)
        super().__init__(L,l, H, h, delta, U, q, Re, N, filestr)

class BFS_H1p5L4_delta0p25(BFS_deltaSmooth):
    def __init__ (self, N):
        L=4
        l=2
        H=1.5
        h=1 
        delta=0.25
        U=0
        q=2 
        Re=0
        name_str = "dBFS_H1p5L4_d0p25_Re0_Q2_U0"
        filestr = "./examples/%s/%s_N%s"%(name_str, name_str, N)
        super().__init__(L, l, H, h, delta, U, q, Re, N, filestr)


class BFS_H1p5L4_delta0p125(BFS_deltaSmooth):
    def __init__ (self, N):
        L=4
        l=2
        H= 1.5 
        h=1 
        delta=0.125
        U=0
        q=2 
        Re=0
        name_str = "dBFS_H1p5L4_d0p125_Re0_Q2_U0"
        filestr = "./examples/%s/%s_N%s"%(name_str, name_str, N)
        super().__init__(L, l, H, h, delta, U, q, Re, N, filestr)


        
class BFS_H1p5L4_delta0(BFS):
    def __init__ (self, N):
        L=4
        l=2
        H= 1.5 
        h=1 

        U=0
        q=2 
        Re=0
        name_str = "dBFS_H1p5L4_d0_Re0_Q2_U0"
        filestr = "./examples/%s/%s_N%s"%(name_str, name_str, N)
        super().__init__(L, l, H, h, U, q, Re, N, filestr)

#------------------------------------------------------------------------------
# <<====== Corner Removal =======>>
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
class BFS_noEddy_Re0_Q2_U0(PWLinear):
    def __init__(self, N):
        x0 = 0
        xf = 4
        y0 = 0
        yf = 2
        y_detatch = 0.4
        x_detatch = 1.35
        x_peaks = [x0,1,x_detatch,xf]
        y_peaks = [[yf,1],[1,y_detatch],[0,0],[0,yf]]
        U = 0
        q = 2
        p0 = 0
        Re = 0
        filestr = "./examples/BFS_H2L4_noEddy_Re0_Q2_U0_N%d"%N
        super().__init__(x0, xf, y0, yf, N, U, q, Re, p0, filestr, x_peaks, y_peaks)
      
#------------------------------------------------------------------------------
# <<====== Half-Hexagon textured slider =======>>
#------------------------------------------------------------------------------


#------------------------------------------------------------------------------
class HexSlider_Re0_Q2_U0(PWLinear):
    def __init__(self, N):
        x0 = 0
        xf = 4
        y0 = 0
        yf = 2
        x_peaks = [x0,1,2,3,xf]
        y_peaks = [[yf,1],[1,0.5],[0,0],[0.5,1],[1,yf]]
        U = 0
        q = 2
        p_amb = 0
        Re = 0
        filestr = "./examples/HexSlider_Re0_Q2_U0_N%d"%N
        super().__init__(x0, xf, y0, yf, N, U, q, Re, p_amb, filestr, x_peaks, y_peaks)


class HexSlider_Re05_Q2_U0(PWLinear):
    def __init__(self, N):
        x0 = 0
        xf = 4
        y0 = 0
        yf = 2
        x_peaks = [x0,1,2,3,xf]
        y_peaks = [[yf,1],[1,0.5],[0,0],[0.5,1],[1,yf]]
        U = 0
        q = 2
        p0 = 0
        Re = 0.5
        filestr = "./examples/HexSlider_Re05_Q2_U0_N%d"%N
        super().__init__(x0, xf, y0, yf, N, U, q, Re, p0, filestr, x_peaks, y_peaks)
        
class HexSlider_Re1_Q2_U0(PWLinear):
    def __init__(self, N):
        x0 = 0
        xf = 4
        y0 = 0
        yf = 2
        x_peaks = [x0,1,2,3,xf]
        y_peaks = [[yf,1],[1,0.5],[0,0],[0.5,1],[1,yf]]
        U = 0
        q = 2
        p0 = 0
        Re = 1
        filestr = "./examples/HexSlider_Re1_Q2_U0_N%d"%N
        super().__init__(x0, xf, y0, yf, N, U, q, Re, p0, filestr, x_peaks, y_peaks)

       
#------------------------------------------------------------------------------
# <<====== Rectangle textured slider =======>>
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
#------------------------------------------------------------------------------
class RectSlider_H2L4_Re0_Q2_U1(RectSlider):
    def __init__(self, N):
        L = 4
        l1 = 2
        l2 = 2
        H = 2
        h = 1
        U = 1
        q = 2
        Re = 0
        filestr = "./examples/slider_rect_Re0_U1_N%d"%N
        super().__init__(L, l1, l2, H, h, U, q, Re, N, filestr)   
        
class RectSlider_H2L4_Re05_Q2_U1(RectSlider):
    def __init__(self, N):
        L = 4
        l1 = 2
        l2 = 2
        H = 2
        h = 1
        U = 1
        q = 2
        Re = 0.5
        filestr = "./examples/slider_rect_Re0_U1_N%d"%N
        super().__init__(L, l1, l2, H, h, U, q, Re, N, filestr)   

#------------------------------------------------------------------------------
# <<====== Triangle slider =======>>
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
class TriSlider_Re05_Q1_U05(PWLinear):
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
        filestr = "./examples/TriSlider_Re0_Q1_U05_N%d"%N
        super().__init__(x0, xf, y0, yf, N, U, q, Re, p0, filestr, x_peaks, y_peaks)
            
    
    
    
    
    
    
    