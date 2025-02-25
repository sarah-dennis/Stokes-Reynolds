 # -*- coding: utf-8 -*-
"""
Created on Tue May 21 16:15:37 2024

@author: sarah
"""

from stokes_heights import PWLinear

#------------------------------------------------------------------------------
# <<=============================== CAVITY ==================================>>
#------------------------------------------------------------------------------
class TriCavity(PWLinear):
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
        name_str = "TriCavity_H4_Re0_Q2_U1"
        filestr = "./examples/%s/%s_N%s"%(name_str, name_str, N)
        super().__init__(x0, xf, y0, yf, N, U, q, Re, filestr, x_peaks, y_peaks)

#------------------------------------------------------------------------------
# <<===================== BACKWARD FACING STEP ==============================>>
#------------------------------------------------------------------------------
class BFS(PWLinear):
    def __init__ (self, L, l, H, h, U, q, Re, N, filestr):
        x0 = 0
        xf = L
        y0 = 0
        yf = H
        x_peaks = [x0, x0+l, xf]
        y_peaks=[[yf,yf-h],[yf-h,0],[0,yf]]
    
        super().__init__(x0, xf, y0, yf, N, U, q, Re, filestr, x_peaks, y_peaks)
        
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
        name_str = "BFS_biswas_Re0"
        filestr = "./examples/%s/%s_N%s"%(name_str, name_str, N)
        super().__init__(L, l, H, h, U, q, Re, N, filestr)

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
        name_str = "BFS_H2p75L4_Re0_Q2_U0"
        filestr = "./examples/%s/%s_N%s"%(name_str, name_str, N)
        super().__init__(L, l, H, h, U, q, Re, N, filestr)        

class BFS_H2p5L4_Re0_Q2_U0(BFS):
    def __init__(self, N):
        L = 4
        l=1
        H = 2.5
        h = 1   
        U = 0
        Re = 0
        q = 2
        name_str = "BFS_H2p5L4_Re0_Q2_U0"
        filestr = "./examples/%s/%s_N%s"%(name_str, name_str, N)
        super().__init__(L, l, H, h, U, q, Re, N, filestr)        

class BFS_H2p25L4_Re0_Q2_U0(BFS):
    def __init__(self, N):
        L = 4
        l=1
        H = 2.25
        h = 1   
        U = 0
        Re = 0
        q = 2
        name_str = "BFS_H2p25L4_Re0_Q2_U0"
        filestr = "./examples/%s/%s_N%s"%(name_str, name_str, N)
        super().__init__(L, l, H, h, U, q, Re, N, filestr)    

class BFS_H2L4_Re0_Q2_U0(BFS):
    def __init__(self, N):
        L = 4
        l = 1
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
class BFS_H2p75L4_Re1_Q2_U0(BFS):
    def __init__(self, N):
        L = 4
        l=1
        H = 2.75
        h = 1   
        U = 0
        Re = 1
        q = 2
        name_str = "BFS_H2p75L4_Re1_Q2_U0"
        filestr = "./examples/%s/%s_N%s"%(name_str, name_str, N)
        super().__init__(L, l, H, h, U, q, Re, N, filestr)
        
class BFS_H2p5L4_Re1_Q2_U0(BFS):
    def __init__(self, N):
        L = 4
        l=1
        H = 2.5
        h = 1   
        U = 0
        Re = 1
        q = 2
        name_str = "BFS_H2p5L4_Re1_Q2_U0"
        filestr = "./examples/%s/%s_N%s"%(name_str, name_str, N)
        super().__init__(L, l, H, h, U, q, Re, N, filestr)

class BFS_H2p25L4_Re1_Q2_U0(BFS):
    def __init__(self, N):
        L = 4
        l=1
        H = 2.25
        h = 1   
        U = 0
        Re = 1
        q = 2
        name_str = "BFS_H2p25L4_Re1_Q2_U0"
        filestr = "./examples/%s/%s_N%s"%(name_str, name_str, N)
        super().__init__(L, l, H, h, U, q, Re, N, filestr)
        
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
        
class BFS_H1p5L4_Re1_Q2_U0(BFS):
    def __init__(self, N):
        L = 4
        l=1
        H = 1.5
        h = 1   
        U = 0
        Re = 1
        q = 2
        name_str = "BFS_H1p5L4_Re1_Q2_U0"
        filestr = "./examples/%s/%s_N%s"%(name_str, name_str, N)
        super().__init__(L, l, H, h, U, q, Re, N, filestr)
  
        
class BFS_H1p25L4_Re1_Q2_U0(BFS):
    def __init__(self, N):
        L = 4
        l=1
        H = 1.25
        h = 1   
        U = 0
        Re = 1
        q = 2
        name_str = "BFS_H1p25L4_Re1_Q2_U0"
        filestr = "./examples/%s/%s_N%s"%(name_str, name_str, N)
        super().__init__(L, l, H, h, U, q, Re, N, filestr)
  
class BFS_H1p125L4_Re1_Q2_U0(BFS):
    def __init__(self, N):
        L = 4
        l=1
        H = 1.125
        h = 1   
        U = 0
        Re = 1
        q = 2
        name_str = "BFS_H1p125L4_Re1_Q2_U0"
        filestr = "./examples/%s/%s_N%s"%(name_str, name_str, N)
        super().__init__(L, l, H, h, U, q, Re, N, filestr)
  
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
        name_str = "BFS_H2p75L4_Re0p5_Q2_U0"
        filestr = "./examples/%s/%s_N%s"%(name_str, name_str, N)
        super().__init__(L, l, H, h, U, q, Re, N, filestr)

class BFS_H2p5L4_Re0p5_Q2_U0(BFS):
    def __init__(self, N):
        L = 4
        l=1
        H = 2.5
        h = 1   
        U = 0
        Re = 0.5
        q = 2
        name_str = "BFS_H2p5L4_Re0p5_Q2_U0"
        filestr = "./examples/%s/%s_N%s"%(name_str, name_str, N)
        super().__init__(L, l, H, h, U, q, Re, N, filestr)
        
class BFS_H2p25L4_Re0p5_Q2_U0(BFS):
    def __init__(self, N):
        L = 4
        l=1
        H = 2.25
        h = 1   
        U = 0
        Re = 0.5
        q = 2
        name_str = "BFS_H2p25L4_Re0p5_Q2_U0"
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
        
class BFS_H1p5L4_Re0p5_Q2_U0(BFS):
    def __init__(self, N):
        L = 4
        l=1
        H = 1.5
        h = 1   
        U = 0
        Re = 0.5
        q = 2
        name_str = "BFS_H1p5L4_Re0p5_Q2_U0"
        filestr = "./examples/%s/%s_N%s"%(name_str, name_str, N)
        super().__init__(L, l, H, h, U, q, Re, N, filestr)
        
class BFS_H1p25L4_Re0p5_Q2_U0(BFS):
    def __init__(self, N):
        L = 4
        l=1
        H = 1.25
        h = 1   
        U = 0
        Re = 0.5
        q = 2
        name_str = "BFS_H1p25L4_Re0p5_Q2_U0"
        filestr = "./examples/%s/%s_N%s"%(name_str, name_str, N)
        super().__init__(L, l, H, h, U, q, Re, N, filestr)
        
class BFS_H1p125L4_Re0p5_Q2_U0(BFS):
    def __init__(self, N):
        L = 4
        l=1
        H = 1.125
        h = 1   
        U = 0
        Re = 0.5
        q = 2
        name_str = "BFS_H1p125L4_Re0p5_Q2_U0"
        filestr = "./examples/%s/%s_N%s"%(name_str, name_str, N)
        super().__init__(L, l, H, h, U, q, Re, N, filestr)
     
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
        name_str = "BFS_H2p75L4_Re0p25_Q2_U0"
        filestr = "./examples/%s/%s_N%s"%(name_str, name_str, N)
        super().__init__(L, l, H, h, U, q, Re, N, filestr)

class BFS_H2p5L4_Re0p25_Q2_U0(BFS):
    def __init__(self, N):
        L = 4
        l=1
        H = 2.5
        h = 1   
        U = 0
        Re = 0.25
        q = 2
        name_str = "BFS_H2p5L4_Re0p25_Q2_U0"
        filestr = "./examples/%s/%s_N%s"%(name_str, name_str, N)
        super().__init__(L, l, H, h, U, q, Re, N, filestr)
        
class BFS_H2p25L4_Re0p25_Q2_U0(BFS):
    def __init__(self, N):
        L = 4
        l=1
        H = 2.25
        h = 1   
        U = 0
        Re = 0.25
        q = 2
        name_str = "BFS_H2p25L4_Re0p25_Q2_U0"
        filestr = "./examples/%s/%s_N%s"%(name_str, name_str, N)
        super().__init__(L, l, H, h, U, q, Re, N, filestr)
       
class BFS_H2L4_Re0p25_Q2_U0(BFS):
    def __init__(self, N):
        L = 4
        l=1
        H = 2
        h = 1   
        U = 0
        Re = 0.25
        q = 2
        name_str = "BFS_H2L4_Re0p25_Q2_U0"
        filestr = "./examples/%s/%s_N%s"%(name_str, name_str, N)
        super().__init__(L, l, H, h, U, q, Re, N, filestr)
        
class BFS_H1p5L4_Re0p25_Q2_U0(BFS):
    def __init__(self, N):
        L = 4
        l=1
        H = 1.5
        h = 1   
        U = 0
        Re = 0.25
        q = 2
        name_str = "BFS_H1p5L4_Re0p25_Q2_U0"
        filestr = "./examples/%s/%s_N%s"%(name_str, name_str, N)
        super().__init__(L, l, H, h, U, q, Re, N, filestr)
        
class BFS_H1p25L4_Re0p25_Q2_U0(BFS):
    def __init__(self, N):
        L = 4
        l=1
        H = 1.25
        h = 1   
        U = 0
        Re = 0.25
        q = 2
        name_str = "BFS_H1p25L4_Re0p25_Q2_U0"
        filestr = "./examples/%s/%s_N%s"%(name_str, name_str, N)
        super().__init__(L, l, H, h, U, q, Re, N, filestr)
        
class BFS_H1p125L4_Re0p25_Q2_U0(BFS):
    def __init__(self, N):
        L = 4
        l=1
        H = 1.125
        h = 1   
        U = 0
        Re = 0.25
        q = 2
        name_str = "BFS_H1p125L4_Re0p25_Q2_U0"
        filestr = "./examples/%s/%s_N%s"%(name_str, name_str, N)
        super().__init__(L, l, H, h, U, q, Re, N, filestr)
         
#------------------------------------------------------------------------------
# <<======================== Smoothed Step  =================================>>
#------------------------------------------------------------------------------
class BFS_deltaSmooth(PWLinear):
    def __init__ (self, L, l,H, h, delta, U, q, Re, N, filestr):
        x0 = 0
        xf = L
        y0 = 0
        yf = H

        x_peaks = [x0, l-delta, l, l+delta, xf]
        y_peaks=[[yf,yf-h],[yf-h,yf-h],[yf-h-(H-h)/2,yf-h-(H-h)/2],[0,0],[0,yf]]

        super().__init__(x0, xf, y0, yf, N, U, q, Re, filestr, x_peaks, y_peaks)

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
        name_str = "dBFS_H2L4_d1p5_Re0_Q2_U0"
        filestr = "./examples/%s/%s_N%s"%(name_str, name_str, N)
        super().__init__(L, l, H, h, delta, U, q, Re, N, filestr)


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
        name_str = "dBFS_H2L4_d1p25_Re0_Q2_U0"
        filestr = "./examples/%s/%s_N%s"%(name_str, name_str, N)
        super().__init__(L, l, H, h, delta, U, q, Re, N, filestr)
        
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
        name_str = "dBFS_H2L4_d1_Re0_Q2_U0"
        filestr = "./examples/%s/%s_N%s"%(name_str, name_str, N)
        super().__init__(L, l, H, h, delta, U, q, Re, N, filestr)

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
        name_str = "dBFS_H2L4_d0p75_Re0_Q2_U0"
        filestr = "./examples/%s/%s_N%s"%(name_str, name_str, N)
        super().__init__(L, l, H, h, delta, U, q, Re, N, filestr)
        

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
        name_str = "dBFS_H2L4_d0p5_Re0_Q2_U0"
        filestr = "./examples/%s/%s_N%s"%(name_str, name_str, N)
        super().__init__(L,l, H, h, delta, U, q, Re, N, filestr)

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
        name_str = "dBFS_H2L4_d0p25_Re0_Q2_U0"
        filestr = "./examples/%s/%s_N%s"%(name_str, name_str, N)
        super().__init__(L, l, H, h, delta, U, q, Re, N, filestr)
        
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
        name_str = "dBFS_H2L4_d0p125_Re0_Q2_U0"
        filestr = "./examples/%s/%s_N%s"%(name_str, name_str, N)
        super().__init__(L, l, H, h, delta, U, q, Re, N, filestr)

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
        name_str = "dBFS_H2L4_d0p01_Re0_Q2_U0"
        filestr = "./examples/%s/%s_N%s"%(name_str, name_str, N)
        super().__init__(L, l, H, h, delta, U, q, Re, N, filestr)


class dBFS_H2L4_d0_Re0_Q2_U0(BFS):
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
        name_str = "dBFS_H1p5L4_d1_Re0_Q2_U0"
        filestr = "./examples/%s/%s_N%s"%(name_str, name_str, N)
        super().__init__(L, l, H, h, delta, U, q, Re, N, filestr)
        

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
        name_str = "dBFS_H1p5L4_d0p75_Re0_Q2_U0"
        filestr = "./examples/%s/%s_N%s"%(name_str, name_str, N)
        super().__init__(L, l, H, h, delta, U, q, Re, N, filestr)
        
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
        name_str = "dBFS_H1p5L4_d0p5_Re0_Q2_U0"
        filestr = "./examples/%s/%s_N%s"%(name_str, name_str, N)
        super().__init__(L,l, H, h, delta, U, q, Re, N, filestr)

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
        name_str = "dBFS_H1p5L4_d0p25_Re0_Q2_U0"
        filestr = "./examples/%s/%s_N%s"%(name_str, name_str, N)
        super().__init__(L, l, H, h, delta, U, q, Re, N, filestr)

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
        name_str = "dBFS_H1p5L4_d0p125_Re0_Q2_U0"
        filestr = "./examples/%s/%s_N%s"%(name_str, name_str, N)
        super().__init__(L, l, H, h, delta, U, q, Re, N, filestr)


class dBFS_H1p5L4_d0_Re0_Q2_U0(BFS):
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
        name_str = "dBFS_H1p25L4_d0p75_Re0_Q2_U0"
        filestr = "./examples/%s/%s_N%s"%(name_str,name_str, N)
        super().__init__(L, l, H, h, delta, U, q, Re, N, filestr)

        
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
        name_str = "dBFS_H1p25L4_d0p5_Re0_Q2_U0"
        filestr = "./examples/%s/%s_N%s"%(name_str,name_str, N)
        super().__init__(L, l, H, h, delta, U, q, Re, N, filestr)

        

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
        name_str = "dBFS_H1p25L4_d0p25_Re0_Q2_U0"
        filestr = "./examples/%s/%s_N%s"%(name_str,name_str, N)
        super().__init__(L, l, H, h, delta, U, q, Re, N, filestr)


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
        name_str = "dBFS_H1p25L4_d0p125_Re0_Q2_U0"
        filestr = "./examples/%s/%s_N%s"%(name_str,name_str, N)
        super().__init__(L, l, H, h, delta, U, q, Re, N, filestr)


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
        name_str = "dBFS_H1p25L4_d0p05_Re0_Q2_U0"
        filestr = "./examples/%s/%s_N%s"%(name_str,name_str, N)
        super().__init__(L, l, H, h, delta, U, q, Re, N, filestr)

class dBFS_H1p25L4_d0_Re0_Q2_U0(BFS):
    def __init__ (self, N):
        L=4
        l=2
        H=1.25 
        h=1 

        U=0
        q=2 
        Re=0
        name_str = "dBFS_H1p25L4_d0_Re0_Q2_U0"
        filestr = "./examples/%s/%s_N%s"%(name_str, name_str, N)
        super().__init__(L, l, H, h, U, q, Re, N, filestr)

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
        p0 = 0
        Re = 0
        name_str = "BFS_H2p75L4_noEddy_Re0_Q2_U0"
        filestr = "./examples/%s/%s_N%s"%(name_str,name_str, N)
        super().__init__(x0, xf, y0, yf, N, U, q, Re, p0, filestr, x_peaks, y_peaks)
        

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
        p0 = 0
        Re = 0
        name_str = "BFS_H2p5L4_noEddy_Re0_Q2_U0"
        filestr = "./examples/%s/%s_N%s"%(name_str,name_str, N)

        super().__init__(x0, xf, y0, yf, N, U, q, Re, p0, filestr, x_peaks, y_peaks)
        

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
        p0 = 0
        Re = 0
        name_str = "BFS_H2p25L4_noEddy_Re0_Q2_U0"
        filestr = "./examples/%s/%s_N%s"%(name_str,name_str, N)

        super().__init__(x0, xf, y0, yf, N, U, q, Re, p0, filestr, x_peaks, y_peaks)
        
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
        p0 = 0
        Re = 0
        name_str = "BFS_H2L4_noEddy_Re0_Q2_U0"
        filestr = "./examples/%s/%s_N%s"%(name_str,name_str, N)

        super().__init__(x0, xf, y0, yf, N, U, q, Re, p0, filestr, x_peaks, y_peaks)


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
        p0 = 0
        Re = 0
        name_str = "BFS_H1p5L4_noEddy_Re0_Q2_U0"
        filestr = "./examples/%s/%s_N%s"%(name_str,name_str, N)

        super().__init__(x0, xf, y0, yf, N, U, q, Re, p0, filestr, x_peaks, y_peaks)      
        
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
        p0 = 0
        Re = 0
        name_str = "BFS_H1p25L4_noEddy_Re0_Q2_U0"
        filestr = "./examples/%s/%s_N%s"%(name_str,name_str, N)

        super().__init__(x0, xf, y0, yf, N, U, q, Re, p0, filestr, x_peaks, y_peaks)      

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
        p0 = 0
        Re = 0
        name_str = "BFS_H2L4_cornerTriA_Re0_Q2_U0"
        filestr = "./examples/%s/%s_N%s"%(name_str,name_str, N)

        super().__init__(x0, xf, y0, yf, N, U, q, Re, p0, filestr, x_peaks, y_peaks)

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
        p0 = 0
        Re = 0
        name_str = "BFS_H2L4_cornerTriB_Re0_Q2_U0"
        filestr = "./examples/%s/%s_N%s"%(name_str,name_str, N)

        super().__init__(x0, xf, y0, yf, N, U, q, Re, p0, filestr, x_peaks, y_peaks)


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
        p0 = 0
        Re = 0
        name_str = "BFS_H2L4_cornerTriC_Re0_Q2_U0"
        filestr = "./examples/%s/%s_N%s"%(name_str,name_str, N)

        super().__init__(x0, xf, y0, yf, N, U, q, Re, p0, filestr, x_peaks, y_peaks)

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
        p0 = 0
        Re = 0
        name_str = "BFS_H2L4_cornerTriD_Re0_Q2_U0"
        filestr = "./examples/%s/%s_N%s"%(name_str,name_str, N)
        super().__init__(x0, xf, y0, yf, N, U, q, Re, p0, filestr, x_peaks, y_peaks)


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
        p0 = 0
        Re = 0
        name_str = "BFS_H2L4_cornerTriE_Re0_Q2_U0"
        filestr = "./examples/%s/%s_N%s"%(name_str,name_str, N)
        super().__init__(x0, xf, y0, yf, N, U, q, Re, p0, filestr, x_peaks, y_peaks)

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
        name_str = 'pipe_Re0_Q2_U0'
        filestr = "./examples/%s/%s_N%d"%(name_str, name_str,N)
        super().__init__(x0, xf, y0, yf, N, U, q, Re, filestr, x_peaks, y_peaks)


    
    
    
    