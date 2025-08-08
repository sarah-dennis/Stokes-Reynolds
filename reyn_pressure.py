 #!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec  5 13:28:42 2024

@author: sarahdennis
"""
import numpy as np


import reyn_pressure_finDiff as fd
from numpy.linalg import solve as np_solve

from reyn_heights import PWL_Height
import reyn_pressure_pwlGMRes as pwlGMRes
from scipy.sparse.linalg import gmres

import reyn_pressure_adjusted

class Pressure:
    def __init__(self, height, U, dP=None, Q=None, ps_1D=None, ps_2D=None):
        # initializing for adjusted solutions generates both ps_1D and ps_2D
    
        # reading in data for any adjusted solution, ps_1D will be missing
        if ps_1D is not None:
            self.ps_1D = ps_1D
        else:
            if dP is not None:
                rhs = fd.make_rhs(height, U, dP)
                mat = fd.make_mat(height)
            elif Q is not None:
                rhs = fd.make_rhs_Q(height, U, Q)
                mat = fd.make_mat_Q(height)
            else:
                raise Exception("Prescribe Q or dP")
            
            self.ps_1D = np_solve(mat, rhs)
        
        # initializing for any 1D Reynolds solution, ps_2D will be missing
        if ps_2D is not None:
            self.ps_2D = ps_2D
        else:
            
            self.ps_2D = self.make_2D_ps(height, ps_1D)
        
        
        self.dP = self.ps_2D[0,-1]- self.ps_2D[0,0]
        
        
    def make_2D_ps(self,height,ps): # sets p(x,y) = p(x) 
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
    def __init__(self, height, U, dP=None, Q=None):
        if dP is not None:
            rhs = fd.make_rhs(height, U, dP)
            mat = fd.make_mat(height)
        elif Q is not None:
            rhs = fd.make_rhs_Q(height, U, Q)
            mat = fd.make_mat_Q(height)
        else:
            raise Exception("Prescribe Q or dP")
        
        ps_1D = np_solve(mat, rhs)
        super().__init__(height, U, dP, Q, ps_1D)


class PwlGMRes_ReynPressure(Pressure):
    def __init__(self, height, U, dP):
            
        if not isinstance(height, PWL_Height):
            raise TypeError('Example is not piecewise linear')
        
        rhs = pwlGMRes.make_rhs(height, U, dP)
        linOp = pwlGMRes.pwlLinOp(height)
        sol_coefs, exit_code = gmres(linOp, rhs, tol=1e-10)
         
        if exit_code != 0:
            raise Exception('gmres did not converge')
            
        ps_1D, Q = pwlGMRes.make_ps(height, sol_coefs, U)
        super().__init__(height, U, dP, Q, ps_1D)
    
class VelAdj_ReynPressure(Pressure):
    def __init__(self, height, U, dP=None, Q=None):
        if dP is not None:
            reyn_pressure = FinDiff_ReynPressure(height, U, dP, Q)
            ps_1D = reyn_pressure.ps_1D
            ps_2D, reyn_derivs, sigma_derivs = reyn_pressure_adjusted.make_adj_ps(height, U, ps_1D, reynFlux=False, TG=False)
        elif Q is not None:
            reyn_pressure = FinDiff_ReynPressure(height, U, dP, Q)
            ps_1D = reyn_pressure.ps_1D
            ps_2D, reyn_derivs, sigma_derivs = reyn_pressure_adjusted.make_adj_ps(height, U, ps_1D, reynFlux=True, TG=False)
        else:
            raise Exception("Prescribe Q or dP")
        
            
        self.reyn_pxs, self.reyn_p2xs, self.reyn_p3xs, self.reyn_p4xs = reyn_derivs
        self.sigmas,self.sigma_xs,self.sigma_2xs = sigma_derivs


        super().__init__(height, U, dP, Q, ps_1D, ps_2D)


class TGAdj_ReynPressure(Pressure):
    def __init__(self, height, U, dP, Q):
        if dP is not None:
            reyn_pressure = FinDiff_ReynPressure(height, U, dP)

            ps_1D = reyn_pressure.ps_1D
            ps_2D, reyn_derivs, _ = reyn_pressure_adjusted.make_adj_ps(height, U, ps_1D, TG=True)
            
        elif Q is not None:
            reyn_pressure = FinDiff_ReynPressure(height, U, Q)

            ps_1D = reyn_pressure.ps_1D
            ps_2D, reyn_derivs, _ = reyn_pressure_adjusted.make_adj_ps(height, U, ps_1D, TG=True)

        self.reyn_pxs, self.reyn_p2xs, self.reyn_p3xs, self.reyn_p4xs = reyn_derivs
       
        super().__init__(height, U, dP, Q, ps_1D, ps_2D)
