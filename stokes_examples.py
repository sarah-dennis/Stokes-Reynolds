 # -*- coding: utf-8 -*-
"""
Created on Tue May 21 16:15:37 2024

@author: sarah
"""
import numpy as np
from stokes_heights import PWLinear
import graphics
#------------------------------------------------------------------------------
# <<=============================== CAVITY ==================================>>
#------------------------------------------------------------------------------
class TriCavity_Re0_U1(PWLinear):
    def __init__ (self, N):
        x0 = 0
        xf = 2
        l=xf//2
        y0 = 0
        yf = 4
        x_peaks = [x0, x0+l, xf]
        y_peaks=[[yf,yf],[0,0],[yf,yf]]
        q=0
        U=1
        Re=0
        namestr = "TriCavity_H4_Re0_Q0_U1"
        super().__init__(x0, xf, y0, yf, N, U, q, Re,namestr, x_peaks, y_peaks)

#------------------------------------------------------------------------------
# <<===================== BACKWARD FACING STEP ==============================>>
#------------------------------------------------------------------------------
class BFS(PWLinear):
    def __init__ (self, L, l, H, h, U, q, Re, N,namestr):
        x0 = 0
        xf = L
        y0 = 0
        yf = H
        x_peaks = [x0, x0+l, xf]
        y_peaks=[[yf,yf-h],[yf-h,0],[0,yf]]
    
        super().__init__(x0, xf, y0, yf, N, U, q, Re,namestr, x_peaks, y_peaks)
        
#------------------------------------------------------------------------------
class BFS_biswas_Re0(BFS):
    def __init__(self, N):
        L = 4
        l = 1
        H = 1.94
        h = 1   
        U = 0
        Re = 1e-4
        q = 1
        namestr = "BFS_biswas_Re0"
        super().__init__(L, l, H, h, U, q, Re, N,namestr)

#------------------------------------------------------------------------------
class BFS_H2p75L4_Re0_Q2_U0(BFS):
    def __init__(self, N):
        L = 4
        l=1
        H = 2.75
        h = 1   
        U = 0
        Re = 0
        q = 2
        namestr = "BFS_H2p75L4_Re0_Q2_U0"
        super().__init__(L, l, H, h, U, q, Re, N,namestr)        

class BFS_H2p5L4_Re0_Q2_U0(BFS):
    def __init__(self, N):
        L = 4
        l=1
        H = 2.5
        h = 1   
        U = 0
        Re = 0
        q = 2
        namestr = "BFS_H2p5L4_Re0_Q2_U0"
        super().__init__(L, l, H, h, U, q, Re, N,namestr)        

class BFS_H2p25L4_Re0_Q2_U0(BFS):
    def __init__(self, N):
        L = 4
        l=1
        H = 2.25
        h = 1   
        U = 0
        Re = 0
        q = 2
        namestr = "BFS_H2p25L4_Re0_Q2_U0"
        super().__init__(L, l, H, h, U, q, Re, N,namestr)    

class BFS_H2L4_Re0_Q2_U0(BFS):
    def __init__(self, N):
        L = 4
        l = 1
        H = 2
        h = 1   
        U = 0
        Re = 0
        q = 2
        namestr = "BFS_H2L4_Re0_Q2_U0"
        super().__init__(L, l, H, h, U, q, Re, N,namestr)
        
class BFS_H2L4_Re0_Q1_U0(BFS):
    def __init__(self, N):
        L = 4
        l = 1
        H = 2
        h = 1   
        U = 0
        Re = 0
        q = 1
        namestr = "BFS_H2L4_Re0_Q1_U0"
        super().__init__(L, l, H, h, U, q, Re, N,namestr)


class BFS_H2L4_Re0_Q1p21_U0(BFS):
    def __init__(self, N):
        L = 4
        l = 1
        H = 2
        h = 1   
        U = 0
        Re = 0
        q = 1.21
        namestr = "BFS_H2L4_Re0_Q1p21_U0"
        super().__init__(L, l, H, h, U, q, Re, N,namestr)


class BFS_H1p5L4_Re0_Q2_U0(BFS):
    def __init__(self, N):
        L = 4
        l=1
        H = 1.5
        h = 1   
        U = 0
        Re = 0
        q = 2
        namestr = "BFS_H1p5L4_Re0_Q2_U0"
        super().__init__(L, l, H, h, U, q, Re, N,namestr)
        
class BFS_H1p25L4_Re0_Q2_U0(BFS):
    def __init__(self, N):
        L = 4
        l = 1
        H = 1.25
        h = 1   
        U = 0
        Re = 0
        q = 2
        namestr = "BFS_H1p25L4_Re0_Q2_U0"
        super().__init__(L, l, H, h, U, q, Re, N,namestr)
         
class BFS_H1p125L4_Re0_Q2_U0(BFS):
    def __init__(self, N):
        L = 4
        l = 1
        H = 1.125
        h = 1   
        U = 0
        Re = 0
        q = 2
        namestr = "BFS_H1p125L4_Re0_Q2_U0"
        super().__init__(L, l, H, h, U, q, Re, N,namestr)

#------------------------------------------------------------------------------
class BFS_H2p75L4_Re1_Q2_U0(BFS):
    def __init__(self, N):
        L = 4
        l=1
        H = 2.75
        h = 1   
        U = 0
        Re = 1
        q = 2
        namestr = "BFS_H2p75L4_Re1_Q2_U0"
        super().__init__(L, l, H, h, U, q, Re, N,namestr)
        
class BFS_H2p5L4_Re1_Q2_U0(BFS):
    def __init__(self, N):
        L = 4
        l=1
        H = 2.5
        h = 1   
        U = 0
        Re = 1
        q = 2
        namestr = "BFS_H2p5L4_Re1_Q2_U0"
        super().__init__(L, l, H, h, U, q, Re, N,namestr)

class BFS_H2p25L4_Re1_Q2_U0(BFS):
    def __init__(self, N):
        L = 4
        l=1
        H = 2.25
        h = 1   
        U = 0
        Re = 1
        q = 2
        namestr = "BFS_H2p25L4_Re1_Q2_U0"
        super().__init__(L, l, H, h, U, q, Re, N,namestr)
        
class BFS_H2L4_Re1_Q2_U0(BFS):
    def __init__(self, N):
        L = 4
        l=1
        H = 2
        h = 1   
        U = 0
        Re = 1
        q = 2
        namestr = "BFS_H2L4_Re1_Q2_U0"
        super().__init__(L, l, H, h, U, q, Re, N,namestr)
        
class BFS_H1p5L4_Re1_Q2_U0(BFS):
    def __init__(self, N):
        L = 4
        l=1
        H = 1.5
        h = 1   
        U = 0
        Re = 1
        q = 2
        namestr = "BFS_H1p5L4_Re1_Q2_U0"
        super().__init__(L, l, H, h, U, q, Re, N,namestr)
  
        
class BFS_H1p25L4_Re1_Q2_U0(BFS):
    def __init__(self, N):
        L = 4
        l=1
        H = 1.25
        h = 1   
        U = 0
        Re = 1
        q = 2
        namestr = "BFS_H1p25L4_Re1_Q2_U0"
        super().__init__(L, l, H, h, U, q, Re, N,namestr)
  
class BFS_H1p125L4_Re1_Q2_U0(BFS):
    def __init__(self, N):
        L = 4
        l=1
        H = 1.125
        h = 1   
        U = 0
        Re = 1
        q = 2
        namestr = "BFS_H1p125L4_Re1_Q2_U0"
        super().__init__(L, l, H, h, U, q, Re, N,namestr)
  
#------------------------------------------------------------------------------
class BFS_H2p75L4_Re0p5_Q2_U0(BFS):
    def __init__(self, N):
        L = 4
        l=1
        H = 2.75
        h = 1   
        U = 0
        Re = 0.5
        q = 2
        namestr = "BFS_H2p75L4_Re0p5_Q2_U0"
        super().__init__(L, l, H, h, U, q, Re, N,namestr)

class BFS_H2p5L4_Re0p5_Q2_U0(BFS):
    def __init__(self, N):
        L = 4
        l=1
        H = 2.5
        h = 1   
        U = 0
        Re = 0.5
        q = 2
        namestr = "BFS_H2p5L4_Re0p5_Q2_U0"
        super().__init__(L, l, H, h, U, q, Re, N,namestr)
        
class BFS_H2p25L4_Re0p5_Q2_U0(BFS):
    def __init__(self, N):
        L = 4
        l=1
        H = 2.25
        h = 1   
        U = 0
        Re = 0.5
        q = 2
        namestr = "BFS_H2p25L4_Re0p5_Q2_U0"
        super().__init__(L, l, H, h, U, q, Re, N,namestr)
       
class BFS_H2L4_Re0p5_Q2_U0(BFS):
    def __init__(self, N):
        L = 4
        l=1
        H = 2
        h = 1   
        U = 0
        Re = 0.5
        q = 2
        namestr = "BFS_H2L4_Re0p5_Q2_U0"
        super().__init__(L, l, H, h, U, q, Re, N,namestr)
        
class BFS_H1p5L4_Re0p5_Q2_U0(BFS):
    def __init__(self, N):
        L = 4
        l=1
        H = 1.5
        h = 1   
        U = 0
        Re = 0.5
        q = 2
        namestr = "BFS_H1p5L4_Re0p5_Q2_U0"
        super().__init__(L, l, H, h, U, q, Re, N,namestr)
        
class BFS_H1p25L4_Re0p5_Q2_U0(BFS):
    def __init__(self, N):
        L = 4
        l=1
        H = 1.25
        h = 1   
        U = 0
        Re = 0.5
        q = 2
        namestr = "BFS_H1p25L4_Re0p5_Q2_U0"
        super().__init__(L, l, H, h, U, q, Re, N,namestr)
        
class BFS_H1p125L4_Re0p5_Q2_U0(BFS):
    def __init__(self, N):
        L = 4
        l=1
        H = 1.125
        h = 1   
        U = 0
        Re = 0.5
        q = 2
        namestr = "BFS_H1p125L4_Re0p5_Q2_U0"
        super().__init__(L, l, H, h, U, q, Re, N,namestr)
     
#------------------------------------------------------------------------------
class BFS_H2p75L4_Re0p25_Q2_U0(BFS):
    def __init__(self, N):
        L = 4
        l=1
        H = 2.75
        h = 1   
        U = 0
        Re = 0.25
        q = 2
        namestr = "BFS_H2p75L4_Re0p25_Q2_U0"
        super().__init__(L, l, H, h, U, q, Re, N,namestr)

class BFS_H2p5L4_Re0p25_Q2_U0(BFS):
    def __init__(self, N):
        L = 4
        l=1
        H = 2.5
        h = 1   
        U = 0
        Re = 0.25
        q = 2
        namestr = "BFS_H2p5L4_Re0p25_Q2_U0"
        super().__init__(L, l, H, h, U, q, Re, N,namestr)
        
class BFS_H2p25L4_Re0p25_Q2_U0(BFS):
    def __init__(self, N):
        L = 4
        l=1
        H = 2.25
        h = 1   
        U = 0
        Re = 0.25
        q = 2
        namestr = "BFS_H2p25L4_Re0p25_Q2_U0"
        super().__init__(L, l, H, h, U, q, Re, N,namestr)
       
class BFS_H2L4_Re0p25_Q2_U0(BFS):
    def __init__(self, N):
        L = 4
        l=1
        H = 2
        h = 1   
        U = 0
        Re = 0.25
        q = 2
        namestr = "BFS_H2L4_Re0p25_Q2_U0"
        super().__init__(L, l, H, h, U, q, Re, N,namestr)
        
class BFS_H1p5L4_Re0p25_Q2_U0(BFS):
    def __init__(self, N):
        L = 4
        l=1
        H = 1.5
        h = 1   
        U = 0
        Re = 0.25
        q = 2
        namestr = "BFS_H1p5L4_Re0p25_Q2_U0"
         
        super().__init__(L, l, H, h, U, q, Re, N,namestr)
        
class BFS_H1p25L4_Re0p25_Q2_U0(BFS):
    def __init__(self, N):
        L = 4
        l=1
        H = 1.25
        h = 1   
        U = 0
        Re = 0.25
        q = 2
        namestr = "BFS_H1p25L4_Re0p25_Q2_U0"
         
        super().__init__(L, l, H, h, U, q, Re, N,namestr)
        
class BFS_H1p125L4_Re0p25_Q2_U0(BFS):
    def __init__(self, N):
        L = 4
        l=1
        H = 1.125
        h = 1   
        U = 0
        Re = 0.25
        q = 2
        namestr = "BFS_H1p125L4_Re0p25_Q2_U0"
         
        super().__init__(L, l, H, h, U, q, Re, N,namestr)
        

class BFS_centered_Re0_Q0p74_U0(BFS):
    def __init__ (self, N):
        L=4
        l=2
        H=2
        h=1
        U=0 
        Re=0 
        q=0.74
        namestr = "cBFS_H1p125L4_Re0p25_Q2_U0"
        super().__init__(L, l, H, h, U, q, Re, N,namestr)
        
        

#------------------------------------------------------------------------------
# <<======================== Smoothed Step  =================================>>
#------------------------------------------------------------------------------
class BFS_deltaSmooth(PWLinear):
    def __init__ (self, L, l,H, h, delta, U, q, Re, N,namestr):
        x0 = 0
        xf = L
        y0 = 0
        yf = H

        x_peaks = [x0, l-delta, l, l+delta, xf]
        y_peaks=[[yf,yf-h],[yf-h,yf-h],[yf-h-(H-h)/2,yf-h-(H-h)/2],[0,0],[0,yf]]

        super().__init__(x0, xf, y0, yf, N, U, q, Re,namestr, x_peaks, y_peaks)

#------------------------------------------------------------------------------

class dBFS_H2L4_d1p5_Re0_Q2_U0(BFS_deltaSmooth):
    def __init__ (self, N):
        L=4
        l=2
        H= 2 
        h=1 
        delta=1.5
        U=0
        q=2
        Re=0
        namestr = "dBFS_H2L4_d1p5_Re0_Q2_U0"
         
        super().__init__(L, l, H, h, delta, U, q, Re, N,namestr)


class dBFS_H2L4_d1p25_Re0_Q2_U0(BFS_deltaSmooth):
    def __init__ (self, N):
        L=4
        l=2
        H=2 
        h=1 
        delta=1.25
        U=0
        q=2
        Re=0
        namestr = "dBFS_H2L4_d1p25_Re0_Q2_U0"
         
        super().__init__(L, l, H, h, delta, U, q, Re, N,namestr)
        
class dBFS_H2L4_d1_Re0_Q2_U0(BFS_deltaSmooth):
    def __init__ (self, N):
        L=4
        l=2
        H= 2 
        h=1 
        delta=1
        U=0
        q=2
        Re=0
        namestr = "dBFS_H2L4_d1_Re0_Q2_U0"
         
        super().__init__(L, l, H, h, delta, U, q, Re, N,namestr)

class dBFS_H2L4_d0p75_Re0_Q2_U0(BFS_deltaSmooth):
    def __init__ (self, N):
        L=4
        l=2
        H= 2 
        h=1 
        delta=0.75
        U=0
        q=2
        Re=0
        namestr = "dBFS_H2L4_d0p75_Re0_Q2_U0"
         
        super().__init__(L, l, H, h, delta, U, q, Re, N,namestr)
        

class dBFS_H2L4_d0p5_Re0_Q2_U0(BFS_deltaSmooth):
    def __init__ (self, N):
        L=4
        l=2
        H= 2 
        h=1 
        delta=0.5
        U=0
        q=2 
        Re=0
        namestr = "dBFS_H2L4_d0p5_Re0_Q2_U0"
         
        super().__init__(L,l, H, h, delta, U, q, Re, N,namestr)

class dBFS_H2L4_d0p25_Re0_Q2_U0(BFS_deltaSmooth):
    def __init__ (self, N):
        L=4
        l=2
        H= 2 
        h=1 
        delta=0.25
        U=0
        q=2 
        Re=0
        namestr = "dBFS_H2L4_d0p25_Re0_Q2_U0"
         
        super().__init__(L, l, H, h, delta, U, q, Re, N,namestr)
        
class dBFS_H2L4_d0p125_Re0_Q2_U0(BFS_deltaSmooth):
    def __init__ (self, N):
        L=4
        l=2
        H= 2 
        h=1 
        delta=0.125
        U=0
        q=2 
        Re=0
        namestr = "dBFS_H2L4_d0p125_Re0_Q2_U0"
         
        super().__init__(L, l, H, h, delta, U, q, Re, N,namestr)

class dBFS_H2L4_d0p01_Re0_Q2_U0(BFS_deltaSmooth):
    def __init__ (self, N):
        L=4
        l=2
        H= 2 
        h=1 
        delta=0.05
        U=0
        q=2 
        Re=0
        namestr = "dBFS_H2L4_d0p01_Re0_Q2_U0"
         
        super().__init__(L, l, H, h, delta, U, q, Re, N,namestr)


class dBFS_H2L4_d0_Re0_Q2_U0(BFS):
    def __init__ (self, N):
        L=4
        l=2
        H= 2 
        h=1 
        U=0
        q=2 
        Re=0
        namestr = "dBFS_H2L4_d0_Re0_Q2_U0"
         
        super().__init__(L, l, H, h, U, q, Re, N,namestr)
    
#------------------------------------------------------------------------------

class dBFS_H1p5L4_d1_Re0_Q2_U0(BFS_deltaSmooth):
    def __init__ (self, N):
        L=4
        l=2
        H=1.5 
        h=1 
        delta=1
        U=0
        q=2
        Re=0
        namestr = "dBFS_H1p5L4_d1_Re0_Q2_U0"
         
        super().__init__(L, l, H, h, delta, U, q, Re, N,namestr)
        

class dBFS_H1p5L4_d0p75_Re0_Q2_U0(BFS_deltaSmooth):
    def __init__ (self, N):
        L=4
        l=2
        H=1.5 
        h=1 
        delta=0.75
        U=0
        q=2
        Re=0
        namestr = "dBFS_H1p5L4_d0p75_Re0_Q2_U0"
         
        super().__init__(L, l, H, h, delta, U, q, Re, N,namestr)
        
class dBFS_H1p5L4_d0p5_Re0_Q2_U0(BFS_deltaSmooth):
    def __init__ (self, N):
        L=4
        l=2
        H=1.5
        h=1 
        delta=0.5
        U=0
        q=2 
        Re=0
        namestr = "dBFS_H1p5L4_d0p5_Re0_Q2_U0"
         
        super().__init__(L,l, H, h, delta, U, q, Re, N,namestr)

class dBFS_H1p5L4_d0p25_Re0_Q2_U0(BFS_deltaSmooth):
    def __init__ (self, N):
        L=4
        l=2
        H=1.5
        h=1 
        delta=0.25
        U=0
        q=2 
        Re=0
        namestr = "dBFS_H1p5L4_d0p25_Re0_Q2_U0"
         
        super().__init__(L, l, H, h, delta, U, q, Re, N,namestr)

class dBFS_H1p5L4_d0p125_Re0_Q2_U0(BFS_deltaSmooth):
    def __init__ (self, N):
        L=4
        l=2
        H= 1.5 
        h=1 
        delta=0.125
        U=0
        q=2 
        Re=0
        namestr = "dBFS_H1p5L4_d0p125_Re0_Q2_U0"
         
        super().__init__(L, l, H, h, delta, U, q, Re, N,namestr)


class dBFS_H1p5L4_d0_Re0_Q2_U0(BFS):
    def __init__ (self, N):
        L=4
        l=2
        H= 1.5 
        h=1 

        U=0
        q=2 
        Re=0
        namestr = "dBFS_H1p5L4_d0_Re0_Q2_U0"
         
        super().__init__(L, l, H, h, U, q, Re, N,namestr)
        
#------------------------------------------------------------------------------
class dBFS_H1p25L4_delta0p75_Re0_Q2_U0(BFS_deltaSmooth):
    def __init__ (self, N):
        L=4
        l=L/2
        H= 1.25
        h=1 
        delta=0.75
        U=0
        q=2
        Re=0
        namestr = "dBFS_H1p25L4_d0p75_Re0_Q2_U0"
         
        super().__init__(L, l, H, h, delta, U, q, Re, N,namestr)

        
class dBFS_H1p25L4_d0p5_Re0_Q2_U0(BFS_deltaSmooth):
    def __init__ (self, N):
        L=4
        l=L/2
        H= 1.25
        h=1 
        delta=0.5
        U=0
        q=2
        Re=0
        namestr = "dBFS_H1p25L4_d0p5_Re0_Q2_U0"
         
        super().__init__(L, l, H, h, delta, U, q, Re, N,namestr)

        

class dBFS_H1p25L4_d0p25_Re0_Q2_U0(BFS_deltaSmooth):
    def __init__ (self, N):
        L=4
        l=L/2
        H= 1.25
        h=1 
        delta=0.25
        U=0
        q=2
        Re=0
        namestr = "dBFS_H1p25L4_d0p25_Re0_Q2_U0"
         
        super().__init__(L, l, H, h, delta, U, q, Re, N,namestr)


class dBFS_H1p25L4_d0p125_Re0_Q2_U0(BFS_deltaSmooth):
    def __init__ (self, N):
        L=4
        l=L/2
        H= 1.25
        h=1 
        delta=0.125
        U=0
        q=2
        Re=0
        namestr = "dBFS_H1p25L4_d0p125_Re0_Q2_U0"
         
        super().__init__(L, l, H, h, delta, U, q, Re, N,namestr)


class dBFS_H1p25L4_d0p05_Re0_Q2_U0(BFS_deltaSmooth):
    def __init__ (self, N):
        L=4
        l=L/2
        H= 1.25
        h=1 
        delta=0.05
        U=0
        q=2
        Re=0
        namestr = "dBFS_H1p25L4_d0p05_Re0_Q2_U0"
         
        super().__init__(L, l, H, h, delta, U, q, Re, N,namestr)

class dBFS_H1p25L4_d0_Re0_Q2_U0(BFS):
    def __init__ (self, N):
        L=4
        l=2
        H=1.25 
        h=1 

        U=0
        q=2 
        Re=0
        namestr = "dBFS_H1p25L4_d0_Re0_Q2_U0"
         
        super().__init__(L, l, H, h, U, q, Re, N,namestr)

#------------------------------------------------------------------------------
# <<=========================  WEDGED CORNER BFS ============================>>
#------------------------------------------------------------------------------

class BFS_H2p75L4_noEddy_Re0_Q2_U0(PWLinear):
    def __init__(self, N):
        x0 = 0
        xf = 4
        y0 = 0
        yf = 2.75
        y_detatch = 0.5
        x_detatch = 1.46
        x_peaks = [x0,1,x_detatch,xf]
        y_peaks = [[yf,yf-1],[yf-1,y_detatch],[0,0],[0,yf]]
        U = 0
        q = 2
          
        Re = 0
        namestr = "BFS_H2p75L4_noEddy_Re0_Q2_U0"
         
        super().__init__(x0, xf, y0, yf, N, U, q, Re ,namestr, x_peaks, y_peaks)
        

class BFS_H2p5L4_noEddy_Re0_Q2_U0(PWLinear):
    def __init__(self, N):
        x0 = 0
        xf = 4
        y0 = 0
        yf = 2.5
        y_detatch = 0.48
        x_detatch = 1.43
        x_peaks = [x0,1,x_detatch,xf]
        y_peaks = [[yf,yf-1],[yf-1,y_detatch],[0,0],[0,yf]]
        U = 0
        q = 2
          
        Re = 0
        namestr = "BFS_H2p5L4_noEddy_Re0_Q2_U0"
         

        super().__init__(x0, xf, y0, yf, N, U, q, Re ,namestr, x_peaks, y_peaks)
        

class BFS_H2p25L4_noEddy_Re0_Q2_U0(PWLinear):
    def __init__(self, N):
        x0 = 0
        xf = 4
        y0 = 0
        yf = 2.25
        y_detatch = 0.44
        x_detatch = 1.34
        x_peaks = [x0,1,x_detatch,xf]
        y_peaks = [[yf,yf-1],[yf-1,y_detatch],[0,0],[0,yf]]
        U = 0
        q = 2
          
        Re = 0
        namestr = "BFS_H2p25L4_noEddy_Re0_Q2_U0"
         

        super().__init__(x0, xf, y0, yf, N, U, q, Re ,namestr, x_peaks, y_peaks)
        
class BFS_H2L4_noEddy_Re0_Q2_U0(PWLinear):
    def __init__(self, N):
        x0 = 0
        xf = 4
        y0 = 0
        yf = 2
        y_detatch = 0.41
        x_detatch = 1.36
        x_peaks = [x0,1,x_detatch,xf]
        y_peaks = [[yf,yf-1],[yf-1,y_detatch],[0,0],[0,yf]]
        U = 0
        q = 2
          
        Re = 0
        namestr = "BFS_H2L4_noEddy_Re0_Q2_U0"
         

        super().__init__(x0, xf, y0, yf, N, U, q, Re ,namestr, x_peaks, y_peaks)


class BFS_H1p5L4_noEddy_Re0_Q2_U0(PWLinear):
    def __init__(self, N):
        x0 = 0
        xf = 4
        y0 = 0
        yf = 1.5
        y_detatch = 0.29
        x_detatch = 1.24
        x_peaks = [x0,1,x_detatch,xf]
        y_peaks = [[yf,yf-1],[yf-1,y_detatch],[0,0],[0,yf]]
        U = 0
        q = 2
          
        Re = 0
        namestr = "BFS_H1p5L4_noEddy_Re0_Q2_U0"
         

        super().__init__(x0, xf, y0, yf, N, U, q, Re ,namestr, x_peaks, y_peaks)      
        
class BFS_H1p25L4_noEddy_Re0_Q2_U0(PWLinear):
    def __init__(self, N):
        x0 = 0
        xf = 4
        y0 = 0
        yf = 1.25
        y_detatch = 0.18
        x_detatch = 1.15
        x_peaks = [x0,1,x_detatch,xf]
        y_peaks = [[yf,yf-1],[yf-1,y_detatch],[0,0],[0,yf]]
        U = 0
        q = 2
          
        Re = 0
        namestr = "BFS_H1p25L4_noEddy_Re0_Q2_U0"
         

        super().__init__(x0, xf, y0, yf, N, U, q, Re ,namestr, x_peaks, y_peaks)      

#------------------------------------------------------------------------------
# <<================= VARIATIONS OF WEDGED CORNER BFS =======================>>
#------------------------------------------------------------------------------

class BFS_H2L4_cornerTriA_Re0_Q2_U0(PWLinear):
    def __init__(self, N):
        x0 = 0
        xf = 4
        y0 = 0
        yf = 2
        
        # y_detatch = 0.41
        # x_detatch = 1.36
        y_detatch = 0.20
        x_detatch = 1.17
        x_peaks = [x0,1,x_detatch,xf]
        y_peaks = [[yf,yf-1],[yf-1,y_detatch],[0,0],[0,yf]]
        U = 0
        q = 2
          
        Re = 0
        namestr = "BFS_H2L4_cornerTriA_Re0_Q2_U0"
         

        super().__init__(x0, xf, y0, yf, N, U, q, Re ,namestr, x_peaks, y_peaks)

class BFS_H2L4_cornerTriB_Re0_Q2_U0(PWLinear):
    def __init__(self, N):
        x0 = 0
        xf = 4
        y0 = 0
        yf = 2
        # y_detatch = 0.41
        # x_detatch = 1.36
        y_detatch = 0.6
        x_detatch = 1.52
        x_peaks = [x0,1,x_detatch,xf]
        y_peaks = [[yf,yf-1],[yf-1,y_detatch],[0,0],[0,yf]]
        U = 0
        q = 2
          
        Re = 0
        namestr = "BFS_H2L4_cornerTriB_Re0_Q2_U0"
         

        super().__init__(x0, xf, y0, yf, N, U, q, Re ,namestr, x_peaks, y_peaks)


class BFS_H2L4_cornerTriC_Re0_Q2_U0(PWLinear):
    def __init__(self, N):
        x0 = 0
        xf = 4
        y0 = 0
        yf = 2
        y_detatch = 0.28
        x_detatch = 1.24
        x_peaks = [x0,1,x_detatch,xf]
        y_peaks = [[yf,yf-1],[yf-1,y_detatch],[0,0],[0,yf]]
        U = 0
        q = 2
          
        Re = 0
        namestr = "BFS_H2L4_cornerTriC_Re0_Q2_U0"
         

        super().__init__(x0, xf, y0, yf, N, U, q, Re ,namestr, x_peaks, y_peaks)

class BFS_H2L4_cornerTriD_Re0_Q2_U0(PWLinear):
    def __init__(self, N):
        x0 = 0
        xf = 4
        y0 = 0
        yf = 2
        y_detatch = 0.28
        x_detatch = 0.24 * 2
        x_peaks = [x0,1,1+x_detatch,xf]
        y_peaks = [[yf,yf-1],[yf-1,y_detatch],[0,0],[0,yf]]
        U = 0
        q = 2
          
        Re = 0
        namestr = "BFS_H2L4_cornerTriD_Re0_Q2_U0"
         
        super().__init__(x0, xf, y0, yf, N, U, q, Re ,namestr, x_peaks, y_peaks)


class BFS_H2L4_cornerTriE_Re0_Q2_U0(PWLinear):
    def __init__(self, N):
        x0 = 0
        xf = 4
        y0 = 0
        yf = 2
        y_detatch = 0.28
        x_detatch = 0.24 * 0.5
        x_peaks = [x0,1,1+x_detatch,xf]
        y_peaks = [[yf,yf-1],[yf-1,y_detatch],[0,0],[0,yf]]
        U = 0
        q = 2
          
        Re = 0
        namestr = "BFS_H2L4_cornerTriE_Re0_Q2_U0"
         
        super().__init__(x0, xf, y0, yf, N, U, q, Re ,namestr, x_peaks, y_peaks)

#------------------------------------------------------------------------------
#basic flat pipe
class pipe_Re0(PWLinear):
    def __init__(self, N):
        x0 = 0
        xf = 4
        y0 = 0
        yf = 1
        x_peaks = [x0,xf]
        y_peaks = [[yf,0],[0,yf]]
        U = 0
        q = 2
        Re = 0
        namestr = 'pipe_Re0_Q2_U0'
         
        super().__init__(x0, xf, y0, yf, N, U, q, Re,namestr, x_peaks, y_peaks)

#------------------------------------------------------------------------------
class bump_Re0(PWLinear):
    def __init__(self, N):
        self.h0 = .1
        self.H=2
        self.slope = 4

        self.L = 2
        self.x0 = -self.L
        self.xf = self.L
        x_peaks = [self.x0 + i/N for i in range (int(1 + 2*self.L*N))]
               
        y_peaks_L = [self.H+self.h0] + [self.h(x) for x in x_peaks[:-1]]
        y_peaks_R =  [self.h(x) for x in x_peaks[1:]] + [self.H+self.h0]

        y_peaks = [[y_L, y_R] for y_L,y_R in np.stack((y_peaks_L,y_peaks_R), axis=1)] 
        y0 = np.min(y_peaks) 
        yf =  np.max(y_peaks)
        # graphics.plot_2D(y_peaks, x_peaks, 'bump', ['x', 'y'])
        U = 1
        q = 1
        Re = 0
        namestr = 'Bump_Re0_Q2_U0'
         
        super().__init__(self.x0, self.xf, y0, yf, N, U, q, Re,namestr, x_peaks, y_peaks)


    def h(self, x):
        return (self.H)*(1-(self.slope/2)*(1+np.cos(np.pi*x/self.L)))

 

class logisticBFS(PWLinear):
    def __init__(self, N, slope, flux, Re, namestr):

        
        self.slope = slope
        
        self.h = 1
        self.H=1
        self.L = 4
        self.center=self.L/2
        self.x0 = 0
        self.xf = self.L
        x_peaks = [self.x0 + i/N for i in range (int(1 + self.L*N))]
               
        y_peaks_L = [self.h+self.H] + [self.h_fun(x) for x in x_peaks[:-1]]
        y_peaks_R =  [self.h_fun(x) for x in x_peaks[1:]] + [self.h+self.H]

        y_peaks = [[y_L, y_R] for y_L,y_R in np.stack((y_peaks_L,y_peaks_R), axis=1)] 
        y0 = np.min(y_peaks) 
        yf =  np.max(y_peaks)
        
        U = 0
        q = flux
         
        super().__init__(self.x0, self.xf, y0, yf, N, U, q, Re,namestr, x_peaks, y_peaks)   

    def h_fun(self, x):
        return self.H - (self.H / ( 1 + np.exp((self.center-x)/self.slope)))
 
class logisticBFS_l1_Re0_Q1p19_U0(logisticBFS):
    def __init__(self, N):
        q=1.19 
        slope=1
        Re=0
        namestr = 'logisticBFS_l1_Re0_Q1p19_U0'
        super().__init__(N, slope, q, Re, namestr)

class logisticBFS_l1_Re1_Q1p19_U0(logisticBFS):
    def __init__(self, N):
        q=1.19 
        slope=1
        Re=1
        namestr = 'logisticBFS_l1_Re1_Q1p19_U0'
        super().__init__(N, slope, q, Re, namestr)
        
class logisticBFS_l0p5_Re0_Q0p99_U0(logisticBFS):
    def __init__(self, N):
        q=0.99 
        slope=1/2
        Re=0
        namestr = 'logisticBFS_l0p5_Re0_Q0p99_U0'
        super().__init__(N, slope, q, Re, namestr)
        
class logisticBFS_l0p5_Re1_Q0p99_U0(logisticBFS):
    def __init__(self, N):
        q=0.99 
        slope=1/2
        Re=1
        namestr = 'logisticBFS_l0p5_Re1_Q0p99_U0'
        super().__init__(N, slope, q, Re, namestr)
        
        
        
class logisticBFS_l0p25_Re0_Q0p85_U0(logisticBFS):
    def __init__(self, N):
        q=0.85
        Re=0
        slope=1/4
        namestr = 'logisticBFS_l0p25_Re0_Q0p85_U0'
        super().__init__(N, slope, q, Re, namestr)
        
class logisticBFS_l0p25_Re1_Q0p85_U0(logisticBFS):
    def __init__(self, N):
        q=0.85
        Re=1
        slope=1/4
        namestr = 'logisticBFS_l0p25_Re1_Q0p85_U0'
        super().__init__(N, slope, q, Re, namestr)

class logisticBFS_l0p125_Re0_Q0p79_U0(logisticBFS):
    def __init__(self, N):
        q=0.79
        Re=0
        slope=1/8
        namestr = 'logisticBFS_l0p125_Re0_Q0p79_U0'
        super().__init__(N, slope, q, Re, namestr)
        
class logisticBFS_l0p125_Re1_Q0p79_U0(logisticBFS):
    def __init__(self, N):
        q=0.79
        Re=1
        slope=1/8
        namestr = 'logisticBFS_l0p125_Re1_Q0p79_U0'
        super().__init__(N, slope, q, Re, namestr)

class logisticBFS_l0p0625_Re0_Q0p77_U0(logisticBFS):
    def __init__(self, N):
        Re=0
        q=0.77
        slope=1/16
        namestr = 'logisticBFS_l0p0625_Re0_Q0p77_U0'
        super().__init__(N, slope, q, Re, namestr)

        

class logisticBFS_l0p03125_Re0_Q0p75_U0(logisticBFS):
    def __init__(self, N):
        q=0.75
        Re=0
        slope=1/32
        namestr = 'logisticBFS_l0p03125_Re0_Q0p75_U0'
        super().__init__(N, slope, q, Re, namestr)

