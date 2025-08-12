# -*- coding: utf-8 -*-
"""
Created on Sun Jul 27 18:10:34 2025

@author: sarah
"""
import numpy as np
from stokes_heights import PWLinear


class BFS(PWLinear):
    def __init__ (self, args, U, Q, Re, N):
        h, H, l, L = args
        x0 = 0
        xf = L
        y0 = 0
        yf = H
        x_peaks = [0, l, L]
        y_peaks=[[yf,yf-h],[yf-h,0],[0,yf]]
        namestr= f'BFS_H{H}L{L}_U{U}_Q{Q}_Re{Re}'
        super().__init__(x0, xf, y0, yf, N, U, Q, Re,namestr, x_peaks, y_peaks)
        
class BFS_pwl(PWLinear):
    def __init__ (self, args, U, Q, Re, N):
        H, L, delta = args
        h=1
        x0 = 0
        xf = L
        y0 = 0
        yf = H
        l = L/2
        x_peaks = [x0, l-delta, l, l+delta, xf]
        y_peaks=[[yf,yf-h],[yf-h,yf-h],[yf-h-(H-h)/2,yf-h-(H-h)/2],[0,0],[0,yf]]
        namestr= f'BFS_pwl_H{H}L{L}d{delta}_U{U}_Q{Q}_Re{Re}'
        super().__init__(x0, xf, y0, yf, N, U, Q, Re,namestr, x_peaks, y_peaks)
     
class Logistic(PWLinear):
    def __init__(self,args, U, Q, Re, N):
        H, h, L, delta = args
        self.delta = delta
        self.h = h
        self.H=H
        self.L = L
        self.center=self.L/2
        self.x0 = 0
        self.xf = self.L
        
        self.h_in = self.h + (self.H-self.h) / ( 1 + np.exp(-self.delta*(self.center-self.x0)))
        
        self.h_out = self.h + (self.H-self.h) / ( 1 + np.exp(-self.delta*(self.center-self.xf)))
        
        self.h_max = max(self.h_in, self.h_out)
        
        x_peaks = [self.x0 + i/N for i in range (int(1 + self.L*N))]
               
        y_peaks_L = [self.h+self.H] + [self.h_fun(x) for x in x_peaks[:-1]]
        y_peaks_R =  [self.h_fun(x) for x in x_peaks[1:]] + [self.h+self.H]

        y_peaks = [[y_L, y_R] for y_L,y_R in np.stack((y_peaks_L,y_peaks_R), axis=1)] 
        y0 = 0
        yf = self.h_max
        
        namestr= f'logistic_H{H}L{L}d{delta}_U{U}_Q{Q}_Re{Re}'
        super().__init__(self.x0, self.xf, y0, yf, N, U, Q, Re,namestr, x_peaks, y_peaks)   

    def h_fun(self, x):
        return self.h_max - (self.h + (self.H-self.h) / ( 1 + np.exp(self.delta*(self.center-x))))
     
class Sinusoid(PWLinear):
    def __init__(self, args, U, Q, Re, N):
        H, h,  L, delta = args
        self.h=h
        self.H=H
        self.k = delta

        self.L = L
        self.x0 = -self.L
        self.xf = self.L
        x_peaks = [self.x0 + i/N for i in range (int(1 + 2*self.L*N))]
        h_peaks = [self.h_fun(x) for x in x_peaks]
        y_peaks_L = [max(h_peaks)+ self.h] + h_peaks[:-1]
        y_peaks_R =  h_peaks[1:] + [max(h_peaks)+ self.h]

        y_peaks = [[y_L, y_R] for y_L,y_R in np.stack((y_peaks_L,y_peaks_R), axis=1)] 
        y0 =  np.min(y_peaks)
        yf =  np.max(y_peaks)

        namestr= f'sinusoid_H{H}h{h}L{L}d{delta}_U{U}_Q{Q}_Re{Re}'
        super().__init__(self.x0, self.xf, y0, yf, N, U, Q, Re,namestr, x_peaks, y_peaks)


    def h_fun(self, x):
        return self.H* (1 + self.h * np.cos(np.pi*self.k*x/self.L))
        
class TriCavity(PWLinear):
    def __init__ (self, args, U, Q, Re, N):
        H, L = args
        x0 = 0
        xf = 2
        l=xf//2
        y0 = 0
        yf = H
        x_peaks = [x0, x0+l, xf]
        y_peaks=[[yf,yf],[0,0],[yf,yf]]
        namestr = f"TriCavity_H{H}L{L}_Re{Re}_Q{Q}_U{U}"
        super().__init__(x0, xf, y0, yf, N, U, Q, Re,namestr, x_peaks, y_peaks)
        
        
class TriSlider(PWLinear):
    def __init__ (self, args, U, Q, Re, N):

        Hin, H,  Hout, Lin, La, Lb, Lout = args

        x0 = 0
        L = Lin+La+Lb+Lout
        xf = L
        y0 = 0
        yf = max(Hin, Hout, H)
        
        # x_peaks = [x0, l-delta, l, l+delta, xf]
        # y_peaks=[[yf,yf-h],[yf-h,yf-h],[yf-h-(H-h)/2,yf-h-(H-h)/2],[0,0],[0,yf]]
        
        x_peaks = [x0, x0+Lin, x0+Lin+La, x0+Lin+La+Lb, xf]
        y_peaks=[[yf,yf-Hin],[yf-Hin,yf-Hin],[yf-H,yf-H],[yf-Hout,yf-Hout], [yf-Hout, yf]]
        namestr = f"TriCavity_H{H}L{L}_Re{Re}_Q{Q}_U{U}"
        super().__init__(x0, xf, y0, yf, N, U, Q, Re,namestr, x_peaks, y_peaks)