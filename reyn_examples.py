# -*- coding: utf-8 -*-
"""
Created on Wed Feb 22 10:01:42 2023

@author: sarah
"""

import numpy as np

from reyn_heights import PWL_Height, SinusoidalHeight, CircleHeight, BumpHeight, LogisticHeight
# -----------------------------------------------------------------------------------------------------------------------------------
#   ____________
#  |_____       |
#        |______|
#


class BFS(PWL_Height):
    def __init__(self, args, N):
        h, H, l, L = args
        x0 = 0
        xf = L
        N_regions = 2
        x_peaks = np.asarray([0, l, L], float)
        h_peaks = np.asarray([[0, h], [h, H], [H, 0]], float)
        namestr = ''#f'BFS_H{int(H)}L{int(xf)}'
        titlestr = ''#f'BFS $H/h={H/h : .3f}$, $L={L : .1f}$'
        super().__init__(x0, xf, N, N_regions, x_peaks, h_peaks, namestr)

# -----------------------------------------------------------------------------------------------------------------------------------
#   ____________
#  |_____       |
#        \______|
#


class BFS_deltaSmooth(PWL_Height):
    def __init__(self, args, N):
        H, delta = args
        l = 2
        h = 1
        x0 = 0
        xf = 4
        L = xf
        N_regions = 4
        x_peaks = np.asarray([x0, l-delta, l, l+delta, xf], float)
        h_peaks = np.asarray(
            [[0, h], [h, h], [h+(H-h)/2, h+(H-h)/2], [H, H], [H, 0]], float)
        namestr = ''#f'dBFS_H{int(H)}L{int(xf)}_d{int(delta)}'
        titlestr = ''#f'$\delta$-BFS $\delta={delta :.3f}$, $H/h={H/h : .3f}$, $L={L : .1f}$'
        super().__init__(x0, xf, N, N_regions, x_peaks, h_peaks, namestr)
# -----------------------------------------------------------------------------------------------------------------------------------
#   ____________
#  |_____       |
#        |      |
#        \______|


class BFS_noEddy(PWL_Height):
    def __init__(self, args, N):
        H = args[0]
        L = 4
        xr = args[1]
        yr = args[2]
        x0 = 0
        xf = L
        x_reattatch = 1 + xr
        y_reattatch = H - yr
        N_regions = 3
        x_peaks = np.asarray([x0, 1, x_reattatch, xf], float)
        h_peaks = np.asarray( [[0, 1], [1, y_reattatch], [H, H], [H, 0]], float)
        namestr = ''#f'cBFS_H{int(H)}L{int(xf)}_xr{int(xr)}yr{int(yr)}'
        titlestr = ''#f'Corner removed BFS $H/h={H/h : .3f}$, $x_r = {x_reattatch:.2f}$, $y_r = {y_reattatch:.2f}$, $L={L : .1f}$ '
        super().__init__(x0, xf, N, N_regions, x_peaks, h_peaks, namestr)

# -----------------------------------------------------------------------------------------------------------------------------------
#   __________
#  |____  ____|
#       \/
#


class TriSlider(PWL_Height):
    def __init__(self, args, N):
        h_in, h, h_out, l_in, l_a, l_b, l_out=args
        N_regions = 4
        
        # x0 = -(l_in + l_a)
        # xf = l_b+l_out
        # x_peaks = np.asarray([x0, -l_a, 0, l_b, xf], float)
        x0 = 0
        xf = x0 + l_in + l_a + l_b + l_out
        x_peaks = np.asarray([x0, x0+l_in, x0+l_in+l_a, x0+l_in+l_a+l_b, xf], float)

        h_peaks = np.asarray(([[0, h_in], [h_in, h_in], [h, h], [h_out, h_out], [h_out, 0]]), float)
        namestr = ''#f'TriSlider_h{int(h)}L{int(xf)}'
        titlestr =''# f'Triangular Textured Slider $h_{in}={h_in:.2f}$, $h_{out}={h_out:.2f}$, $h={h:.2f}$, $l_{a}={l_a:.2f}$, $l_b={l_b:.2f}$'
        super().__init__(x0, xf, N, N_regions, x_peaks, h_peaks, namestr)


# -----------------------------------------------------------------------------------------------------------------------------------
#   ______
#   \    /
#    \  /
#     \/


class TriCavity(PWL_Height):
    def __init__(self, args, N):
        H, l_a, l_b = args
        x0 = 0
        xf = x0 + l_a + l_b
        N_regions = 2

        H = args[0]

        x_peaks = np.asarray([x0, x0 + l_a, xf], float)

        dx = 0.01
        # h = max(args[2],dx) #H/2  # dx
        
        h_peaks = np.asarray(([[0, dx], [H, H], [dx, 0]]), float)
        namestr ='' #'TriCavity_H{int(H)}L{int(xf)}'
        titlestr = 'Triangular Cavity'
        super().__init__(x0, xf, N, N_regions, x_peaks, h_peaks, namestr)
        
# -----------------------------------------------------------------------------------------------------------------------------------
# NOT PIECEWISE LINEAR EXAMPLES...
#-----------------------------------------------------------------------------------------------------------------------------------

class Cylinder(CircleHeight):
    def __init__(self, args, N):
        r = args[0]
        h0 = args[1]
        l = args[2]
        x0 = -(l+r)
        drdx = args[3]
        xf = l+r
        namestr = ''
        # namestr = f'Circ_r{int(r)}h{int(h0)}L{int(xf)}'
        titlestr= 'Cylinder'
        super().__init__(x0, xf, N, r,drdx, h0, namestr)
        
#-----------------------------------------------------------------------------------------------------------------------------------

class Sinusoid(SinusoidalHeight):
    def __init__(self, args, N):
        x0 = 0
        xf = args[2]

        H = args[0]
        h = args[1]
        namestr = ''
        # namestr = f'Sinusoid_H{H}h{h}'
        super().__init__(x0, xf, N, H, h, namestr)

#-----------------------------------------------------------------------------------------------------------------------------------
class LambdaBump(BumpHeight):
    def __init__(self, args, N):
        x0 = -args[2]
        xf = args[2]

        H = args[1]
        lam = args[0]
        h0 = args[3]
        
        namestr = ''#f'Bump_lambda{int(lam)}H{int(H)}'
        super().__init__(x0, xf, N, lam, H, h0, namestr)

#-----------------------------------------------------------------------------------------------------------------------------------
class Logistic(LogisticHeight):
    def __init__(self, args, N):
        H_in, H_out, L, delta = args

        # x0 = -L//2
        # xf = L//2
        # center = 0
        x0 = 0
        xf = L
        center = L//2
        namestr = '' #f'Logistic{delta:.2f}H{H:.2f}'
        titlestr =''# f'Logistic Step'#' $\lambda={delta}$, $dh/dx={-delta/4}$'
        super().__init__(x0, xf, N, H_in, H_out, center, delta, namestr)


