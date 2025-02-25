# -*- coding: utf-8 -*-
"""
Created on Wed Feb 22 10:01:42 2023

@author: sarah
"""

import numpy as np

from reyn_heights import PWL_Height, SinusoidalHeight, CircleHeight

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
        name_str = f'BFS_H{H}L{xf}_dP{dP}_U{U}'
        filestr = f'./reyn_examples/{name_str}/{name_str}_N{N}'
        super().__init__(x0, xf, N, N_regions, x_peaks, h_peaks, U, dP, filestr)
        
class BFS_deltaSmooth(PWL_Height):
    def __init__ (self, U, dP, N, args):
        H = args[0]
        delta = args[1]

        l=2
        h=1
        x0 = 0
        xf = 4
        N_regions = 4
        x_peaks = np.asarray([x0, l-delta, l, l+delta, xf], float)
        h_peaks=np.asarray([[h,h],[h,h],[h+(H-h)/2,h+(H-h)/2],[H,H],[H,H]], float)
        name_str = f'dBFS_H{H}L{xf}_d{delta}_dP{dP}_U{U}'
        filestr = f'./reyn_examples/{name_str}/{name_str}_N{N}'
        super().__init__(x0, xf, N, N_regions, x_peaks, h_peaks, U, dP, filestr)

class BFS_noEddy(PWL_Height):
    def __init__(self, U, dP,N, args):
        H=args[0]
        L=4
        xr=args[1]
        yr=args[2]
        x0 = 0
        xf = L
        x_reattatch=1 +xr
        y_reattatch=H -yr
        N_regions = 3
        x_peaks = np.asarray([x0, 1, x_reattatch, xf], float)
        h_peaks = np.asarray([[1,1],[1,y_reattatch],[y_reattatch,H],[H,H]], float)
        name_str = f'cBFS_H{H}L{xf}_xr{xr}yr{yr}_dP{dP}_U{U}'
        filestr = f'./reyn_examples/{name_str}/{name_str}_N{N}'
        super().__init__(x0, xf, N, N_regions, x_peaks, h_peaks, U, dP, filestr)
        
#-----------------------------------------------------------------------------------------------------------------------------------
# args=[r,k]
class Bump(SinusoidalHeight):
    def __init__(self, U, dP, N, args):
        x0=0
        xf=4
        h_avg=1
        
        r=args[0]
        k=args[1]
        name_str = f'Bump_r{r}k{k}_dP{dP}_U{U}'
        filestr = f'./reyn_examples/{name_str}/{name_str}_N{N}'
        super().__init__(x0, xf, N, h_avg, r, k, U, dP, filestr)

class CircCavity(CircleHeight):
    def __init__(self, U, dP, N, args):
        x0 = 0
        r = args[0]
        h0=args[1]
        l=args[2]
        xf = x0 + 2*r + 2*l
        name_str = f'Circ_r{r}h{h0}L{xf}_dP{dP}_U{U}'
        filestr = f'./reyn_examples/{name_str}/{name_str}_N{N}'
        super().__init__(x0, xf, N, r, h0, l, U, dP, filestr)
    
#-----------------------------------------------------------------------------------------------------------------------------------

class HexSlider(PWL_Height):
    def __init__(self, U, dP,N, args=None):
        x0 = 0
        xf = 4
        N_regions = 4
        x_peaks = np.asarray([x0, 1, 2, 3, xf], float)
        h_peaks = np.asarray(([[0,1],[1,1.5],[2,2],[1.5,1],[1,0]]), float)
        name_str = f'Hex_dP{dP}_U{U}'
        filestr = f'./reyn_examples/{name_str}/{name_str}_N{N}'
        super().__init__(x0, xf, N, N_regions, x_peaks, h_peaks, U, dP, filestr)
        
#-----------------------------------------------------------------------------------------------------------------------------------        

class TriSlider(PWL_Height):
    def __init__(self, U, dP,N, args=None):
        x0 = 0
        xf = 4
        N_regions = 4
        
        x_peaks = np.asarray([x0, 1, 2, 3, xf], float)
        
        h_peaks = np.asarray(([[0,1],[1,1],[2,2],[1,1],[1,0]]), float)
        name_str = f'Tri_dP{dP}_U{U}'
        filestr = f'./reyn_examples/{name_str}/{name_str}_N{N}'
        super().__init__(x0, xf, N, N_regions, x_peaks, h_peaks, U, dP, filestr)
    

class TriCavity(PWL_Height):
    def __init__(self, U, dP,N, args):
        x0 = 0
        xf = x0 + 2*args[1]
        N_regions = 2

        H = args[0]
        x_peaks = np.asarray([x0, xf//2, xf], float)
        
        dx = 1/((xf-x0)*N)
        h_peaks = np.asarray(([[0,0+dx/10],[H,H],[0+dx/10,0]]), float)
        name_str = f'TriCavity_H{int(H)}L{int(xf)}_U{int(U)}'
        filestr = f'./reyn_examples/{name_str}/{name_str}_N{N}'
        super().__init__(x0, xf, N, N_regions, x_peaks, h_peaks, U, dP, filestr)
    












