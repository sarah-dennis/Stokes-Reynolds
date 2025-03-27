# -*- coding: utf-8 -*-
"""
Created on Wed Feb 22 10:01:42 2023

@author: sarah
"""

import numpy as np

from reyn_heights import PWL_Height, SinusoidalHeight, CircleHeight, BumpHeight
#-----------------------------------------------------------------------------------------------------------------------------------
#   ____________
#  |_____       | 
#        |______|
#
class BFS(PWL_Height):
    def __init__(self, U, dP, N, args):
        H = args[0]
        l = args[1]
        h = 1
        x0 = 0
        xf = args[2]
        N_regions = 2
        x_peaks = np.asarray([x0, l, xf],float)
        h_peaks = np.asarray([[0,h],[h,H],[H,0]],float)
        namestr = f'BFS_H{int(H)}L{int(xf)}_dP{int(dP)}_U{int(U)}'
        super().__init__(x0, xf, N, N_regions, x_peaks, h_peaks, U, dP, namestr)

#-----------------------------------------------------------------------------------------------------------------------------------
#   ____________
#  |_____       |
#        \______|
#

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
        namestr = f'dBFS_H{int(H)}L{int(xf)}_d{int(delta)}_dP{int(dP)}_U{int(U)}'
        
        super().__init__(x0, xf, N, N_regions, x_peaks, h_peaks, U, dP, namestr)
#-----------------------------------------------------------------------------------------------------------------------------------
#   ____________
#  |_____       |
#        |      |
#        \______|

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
        namestr = f'cBFS_H{int(H)}L{int(xf)}_xr{int(xr)}yr{int(yr)}_dP{int(dP)}_U{int(U)}'
        super().__init__(x0, xf, N, N_regions, x_peaks, h_peaks, U, dP, namestr)
        

#-----------------------------------------------------------------------------------------------------------------------------------
#   ______________ 
#  |_____    _____|
#        \__/
#
class HexSlider(PWL_Height):
    def __init__(self, U, dP,N, args=None):
        x0 = 0
        xf = 4
        N_regions = 4
        x_peaks = np.asarray([x0, 1, 2, 3, xf], float)
        h_peaks = np.asarray(([[0,1],[1,1.5],[2,2],[1.5,1],[1,0]]), float)
        namestr = f'Hex_dP{int(dP)}_U{int(U)}'
        super().__init__(x0, xf, N, N_regions, x_peaks, h_peaks, U, dP, namestr)
        
#-----------------------------------------------------------------------------------------------------------------------------------        
#   __________ 
#  |____  ____|
#       \/
#
class TriSlider(PWL_Height):
    def __init__(self, U, dP,N, args=None):
        x0 = 0
        xf = 4
        N_regions = 4
        
        x_peaks = np.asarray([x0, 1, 2, 3, xf], float)
        
        h_peaks = np.asarray(([[0,1],[1,1],[2,2],[1,1],[1,0]]), float)
        namestr = f'Tri_dP{int(dP)}_U{int(U)}'
        super().__init__(x0, xf, N, N_regions, x_peaks, h_peaks, U, dP, namestr)

#-----------------------------------------------------------------------------------------------------------------------------------        
#   ______
#   \    /
#    \  /
#     \/

class TriCavity(PWL_Height):
    def __init__(self, U, dP,N, args):
        x0 = 0
        xf = x0 + 2*args[1]
        N_regions = 2 

        H = args[0]
        x_peaks = np.asarray([x0, xf//2, xf], float)
        
        dx = 1/((xf-x0)*N)
        h_peaks = np.asarray(([[0,0+dx],[H,H],[0+dx,0]]), float)
        namestr = f'TriCavity_H{int(H)}L{int(xf)}_U{int(U)}'
        super().__init__(x0, xf, N, N_regions, x_peaks, h_peaks, U, dP, namestr)
    
#-----------------------------------------------------------------------------------------------------------------------------------

class Sinusoid(SinusoidalHeight):
    def __init__(self, U, dP, N, args):
        x0=0
        xf=args[2]

        
        r=args[0]
        h_avg=2*r + r/2
        k=args[1]
        namestr = f'Sinusoid_r{int(r)}k{int(k)}_dP{int(dP)}_U{int(U)}'
        super().__init__(x0, xf, N, h_avg, r, k, U, dP, namestr)

class LambdaBump(BumpHeight):
    def __init__(self, U, dP, N, args):
        x0=-args[2]
        xf=args[2]

        H = args[1]
        lam = args[0]
        namestr = f'Bump_lambda{int(lam)}H{int(H)}_dP{int(dP)}_U{int(U)}'
        super().__init__(x0, xf, N, lam, H, U, dP, namestr)


class Cylinder(CircleHeight):
    def __init__(self, U, dP, N, args):
        
        r = args[0]
        h0 = args[1]
        l=args[2]
        x0 = -(l+r)

        xf = l+r
        namestr = f'Circ_r{int(r)}h{int(h0)}L{int(xf)}_dP{int(dP)}_U{int(U)}'
        super().__init__(x0, xf, N, r, h0, U, dP, namestr)
    









