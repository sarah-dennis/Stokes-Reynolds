# -*- coding: utf-8 -*-
"""
Created on Wed Feb 22 10:01:42 2023

@author: sarah
"""

import numpy as np

from reyn_heights import PWL_Height, SinusoidalHeight

class BFS(PWL_Height):
    def __init__(self, U, dP, N, args):
        H = args[0]
        l = args[1]
        h = 1
        x0 = 0
        xf = 4
        N_regions = 2
        x_peaks = np.asarray([x0, l, xf],float)
        h_peaks = np.asarray([[0,h],[h,H],[H,0]],float)
        super().__init__(x0, xf, N, N_regions, x_peaks, h_peaks,U,dP)
        
class BFS_deltaSmooth(PWL_Height):
    def __init__ (self, U, dP, N, args):
        H = args[0]
        delta = args[1]

        l=2
        h=1
        x0 = 0
        xf = 4
        N_regions = 4
        x_peaks = np.asarray([x0, l-delta, l, l+delta, xf],float)
        h_peaks=np.asarray([[h,h],[h,h],[h+(H-h)/2,h+(H-h)/2],[H,H],[H,H]],float)
        super().__init__(x0, xf, N, N_regions,x_peaks, h_peaks,U,dP)

class BFS_noEddy(PWL_Height):
    def __init__(self, U, dP,N, args):
        H=args[0]
        L=4
        xL=args[1]
        yL=args[2]
        x0 = 0
        xf = L
        x_reattatch=1 +xL
        y_reattatch=H -yL
        N_regions = 3
        x_peaks = np.asarray([x0, 1, x_reattatch, xf],float)
        h_peaks = np.asarray([[1,1],[1,y_reattatch],[y_reattatch,H],[H,H]],float)
        super().__init__(x0, xf, N, N_regions, x_peaks, h_peaks,U,dP)
        
#-----------------------------------------------------------------------------------------------------------------------------------
# args=[r,k]
class Bump(SinusoidalHeight):
    def __init__(self, U, dP, N, args):
        x0=0
        xf=4
        h_avg=1
        
        r=args[0]
        k=args[1]
        super().__init__(x0, xf, N, h_avg, r, k, U, dP)

#-----------------------------------------------------------------------------------------------------------------------------------

class HexSlider(PWL_Height):
    def __init__(self, U, dP,N, args=None):
        x0 = 0
        xf = 4
        N_regions = 4
        x_peaks = np.asarray([x0, 1, 2, 3, xf],float)
        h_peaks = np.asarray(([[0,1],[1,1.5],[2,2],[1.5,1],[1,0]]),float)
        super().__init__(x0, xf, N, N_regions, x_peaks, h_peaks,U,dP)
        
#-----------------------------------------------------------------------------------------------------------------------------------        

class TriSlider(PWL_Height):
    def __init__(self, U, dP,N, args=None):
        x0 = 0
        xf = 4
        N_regions = 4
        
        x_peaks = np.asarray([x0, 1, 2, 3, xf],float)
        
        h_peaks = np.asarray(([[0,1],[1,1],[2,2],[1,1],[1,0]]),float)
        
        super().__init__(x0, xf, N, N_regions, x_peaks, h_peaks,U,dP)
    

class TriCavity(PWL_Height):
    def __init__(self, U, dP,N, args):
        x0 = 0
        xf = 2
        N_regions = 2
        H_mult=args[0]
        H = H_mult*(xf-x0)
        x_peaks = np.asarray([x0, xf//2, xf],float)
        
        dx = 1/((xf-x0)*N)
        h_peaks = np.asarray(([[0,0+dx/10],[H,H],[0+dx/10,0]]),float)
        
        super().__init__(x0, xf, N, N_regions, x_peaks, h_peaks,U,dP)
    
         
class variableSlider(PWL_Height):
    def __init__(self, U, dP,N, args=None):
        x0 = 0
        xf = 8
        N_regions = 6
        
        x_peaks = np.asarray([x0, 1.5, 2.5, 4, 5.5, 7, xf],float)
        
        h_peaks = np.asarray(([[1,1],[1,1],[2,2],[1.5,1.5],[2,2],[1,1],[1,1]]),float)
        
        super().__init__(x0, xf, N, N_regions, x_peaks, h_peaks,U,dP)















