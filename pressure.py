# -*- coding: utf-8 -*-
"""
Created on Mon Jan 26 20:20:09 2026

@author: sarah
"""
import numpy as np
import reyn_pressure_ELT as p_ELT

class Pressure:
    def __init__(self, ps_1D, ps_2D):
        # initializing for adjusted solutions generates both ps_1D and ps_2D
    
        self.ps_1D = ps_1D
        
        self.ps_2D = ps_2D 
    
    def get_dP(self, height):
        if self.ps_2D is None:
            dP =  self.ps_1D[0] - self.ps_1D[-1]

        else:
            ps_2D = np.nan_to_num(self.ps_2D)
            dP = (sum(ps_2D[:,0])/height.hs[0] - sum(ps_2D[:,-1])/height.hs[-1])*height.dy
        return dP

class Reyn_Pressure(Pressure):
    def __init__(self, height, ps_1D):
        ps_2D = self.make_2D_ps(height, ps_1D)
        
        super().__init__(ps_1D, ps_2D)
        
    def make_2D_ps(self, height, ps_1D): # p(x,y) = p(x) 
         ps_2D = np.zeros((height.Ny, height.Nx))
         
         for i in range(height.Nx):
             for j in range(height.Ny):
                 
                 y = height.ys[j]
                 if y <= height.hs[i]:
                     ps_2D[j,i] = ps_1D[i]
                 else:
                     ps_2D[j,i] = None
                
         return ps_2D     

class VA_ELT_Pressure(Pressure):
    def __init__(self, height, BC, reyn_pressure):
        
        ps_1D = reyn_pressure.ps_1D
        ps_2D, reyn_derivs, sigma_derivs = p_ELT.make_ELT_ps(height, BC, ps_1D, TG=False)
    
        self.reyn_pxs, self.reyn_p2xs, self.reyn_p3xs, self.reyn_p4xs = reyn_derivs
        self.sigmas,self.sigma_xs,self.sigma_2xs = sigma_derivs


        super().__init__(ps_1D, ps_2D)


class TG_ELT_Pressure(Pressure):
    def __init__(self, height, BC, reyn_pressure):

        ps_1D = reyn_pressure.ps_1D
        ps_2D, reyn_derivs, sigma_derivs = p_ELT.make_ELT_ps(height, BC, ps_1D, TG=True)
    

        self.reyn_pxs, self.reyn_p2xs, self.reyn_p3xs, self.reyn_p4xs = reyn_derivs
        self.sigmas, self.sigma_xs, self.sigma_2xs = sigma_derivs
        super().__init__(ps_1D, ps_2D)
