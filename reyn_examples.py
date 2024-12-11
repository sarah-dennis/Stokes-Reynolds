# -*- coding: utf-8 -*-
"""
Created on Wed Feb 22 10:01:42 2023

@author: sarah
"""

import numpy as np

from reyn_heights import PWL_Height

class BFS(PWL_Height):
    def __init__(self, U, dP, N):
        H=1.25
        L=4
        h=1
        l=2
        x0 = 0
        xf = L
        N_regions = 2
        x_peaks = np.asarray([x0, l, xf],float)
        h_peaks = np.asarray([[0,h],[h,H],[H,0]],float)
        super().__init__(x0, xf, N, N_regions, x_peaks, h_peaks,U,dP)
        
class BFS_deltaSmooth(PWL_Height):
    def __init__ (self, U, dP, N):
        x0 = 0
        xf = 4
        l=2
        H=1.25
        h=1
        delta = 0.125
        N_regions = 4
        x_peaks = np.asarray([x0, l-delta, l, l+delta, xf],float)
        h_peaks=np.asarray([[h,h],[h,h],[h+(H-h)/2,h+(H-h)/2],[H,H],[H,H]],float)
        super().__init__(x0, xf, N, N_regions,x_peaks, h_peaks,U,dP)

class BFS_noEddy(PWL_Height):
    def __init__(self, U, dP,N):
        H=1.25
        L=4
        xL=0.15
        yL=0.17
        x0 = 0
        xf = L
        x_reattatch=1 +xL
        y_reattatch=H -yL
        N_regions = 3
        x_peaks = np.asarray([x0, 1, x_reattatch, xf],float)
        h_peaks = np.asarray([[1,1],[1,y_reattatch],[y_reattatch,H],[H,H]],float)
        super().__init__(x0, xf, N, N_regions, x_peaks, h_peaks,U,dP)
        
        

#-----------------------------------------------------------------------------------------------------------------------------------

class HexSlider(PWL_Height):
    def __init__(self, U, dP,N):
        x0 = 0
        xf = 4
        N_regions = 4
        x_peaks = np.asarray([x0, 1, 2, 3, xf],float)
        h_peaks = np.asarray(([[1,1],[1,1.5],[2,2],[1.5,1],[1,1]]),float)
        super().__init__(x0, xf, N, N_regions, x_peaks, h_peaks,U,dP)
        
#-----------------------------------------------------------------------------------------------------------------------------------        

class TriSlider(PWL_Height):
    def __init__(self, U, dP,N):
        x0 = 0
        xf = 4
        N_regions = 4
        
        x_peaks = np.asarray([x0, 1, 2, 3, xf],float)
        
        h_peaks = np.asarray(([[1,1],[1,1],[2,2],[1,1],[1,1]]),float)
        
        super().__init__(x0, xf, N, N_regions, x_peaks, h_peaks,U,dP)
        
         
class variableSlider(PWL_Height):
    def __init__(self, U, dP,N):
        x0 = 0
        xf = 8
        N_regions = 6
        
        x_peaks = np.asarray([x0, 1.5, 2.5, 4, 5.5, 7, xf],float)
        
        h_peaks = np.asarray(([[1,1],[1,1],[2,2],[1.5,1.5],[2,2],[1,1],[1,1]]),float)
        
        super().__init__(x0, xf, N, N_regions, x_peaks, h_peaks,U,dP)















