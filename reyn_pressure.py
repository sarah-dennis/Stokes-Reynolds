 #!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec  5 13:28:42 2024

@author: sarahdennis
"""
import numpy as np



from numpy.linalg import solve as np_solve


import reyn_pressure_finDiff as fd
import reyn_pressure_pwlGMRes as pwlGMRes
import reyn_pressure_schur as pwcSchur
import reyn_pressure_adjusted


class Pressure:
    def __init__(self, height, BC, ps_1D=None, ps_2D=None):
        # initializing for adjusted solutions generates both ps_1D and ps_2D
    
        if ps_1D is not None:
            self.ps_1D = ps_1D
        else:
            rhs = fd.make_rhs(height, BC)
            mat = fd.make_mat(height, BC)

            self.ps_1D = np_solve(mat, rhs)
        
        # initializing for any 1D Reynolds solution, ps_2D will be missing
        if ps_2D is not None:
            self.ps_2D = ps_2D
        else:
            self.ps_2D = self.make_2D_ps(height, ps_1D)
        
        
        # self.dP = self.ps_2D[0,-1]- self.ps_2D[0,0]
        self.dP = self.ps_1D[-1]-self.ps_1D[0]
        
        
    def make_2D_ps(self,height,ps): # p(x,y) = p(x) 
        ps_2D = np.zeros((height.Ny, height.Nx))
        
        for i in range(height.Nx):
            for j in range(height.Ny):
                
                y = height.ys[j]
                if y <= height.hs[i]:
                    ps_2D[j,i] = ps[i]
                else:
                    ps_2D[j,i] = None
        return ps_2D
                    
class FinDiff_ReynPressure(Pressure):
    def __init__(self, height, BC):
        ps_1D = fd.fd_solve(height, BC)
        
        super().__init__(height, BC, ps_1D=ps_1D)


class PwlGMRes_ReynPressure(Pressure):
    def __init__(self, height, BC):
            

            
        ps_1D = pwlGMRes.gmres_solve(height, BC)
        super().__init__(height, BC, ps_1D)
    
class Schur_ReynPressure(Pressure):
    def __init__(self, height, BC): 
        
        # ps_1D = pwcSchur.schur_solve_parallel(height, BC)

        ps_1D = pwcSchur.schur_solve(height, BC)

        super().__init__(height, BC, ps_1D)
        
class Schur_parallel_ReynPressure(Pressure):
    def __init__(self, height, BC): 
        
        ps_1D = pwcSchur.schur_solve_parallel(height, BC)

        # ps_1D = pwcSchur.schur_solve(height, BC)

        super().__init__(height, BC, ps_1D)

class VelAdj_ReynPressure(Pressure):
    def __init__(self, height, BC):
        
        reyn_pressure = FinDiff_ReynPressure(height, BC)
        ps_1D = reyn_pressure.ps_1D
        ps_2D, reyn_derivs, sigma_derivs = reyn_pressure_adjusted.make_adj_ps(height, BC, ps_1D, TG=False)
    
        self.reyn_pxs, self.reyn_p2xs, self.reyn_p3xs, self.reyn_p4xs = reyn_derivs
        self.sigmas,self.sigma_xs,self.sigma_2xs = sigma_derivs


        super().__init__(height, BC, ps_1D, ps_2D)


class TGAdj_ReynPressure(Pressure):
    def __init__(self, height, BC):

        reyn_pressure = FinDiff_ReynPressure(height, BC)

        ps_1D = reyn_pressure.ps_1D
        ps_2D, reyn_derivs, _ = reyn_pressure_adjusted.make_adj_ps(height, BC, ps_1D, TG=True)
    

        self.reyn_pxs, self.reyn_p2xs, self.reyn_p3xs, self.reyn_p4xs = reyn_derivs
       
        super().__init__(height, BC, ps_1D, ps_2D)
