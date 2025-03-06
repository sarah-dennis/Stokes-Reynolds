 #!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec  5 13:28:42 2024

@author: sarahdennis
"""
import numpy as np
import reyn_pressure_adjusted

import reyn_pressures_finDiff as fd
from numpy.linalg import solve as np_fd_solve
from scipy.sparse.linalg import gmres as sp_pwl_gmres
from reyn_heights import PWL_Height
import reyn_pressures_pwlinear as pwl

class Pressure:
    def __init__(self, height, ps_1D=None, ps_2D=None):
        # initializing for adjusted solutions generates both ps_1D and ps_2D
    
        # reading in data for any adjusted solution, ps_1D will be missing
        if ps_1D is None:
            rhs = fd.make_rhs(height)
            mat = fd.make_mat(height)
            self.ps_1D = np_fd_solve(mat, rhs)
        else:
            self.ps_1D = ps_1D
        
        # initializing for any 1D Reynolds solution, ps_2D will be missing
        if ps_2D is None:
            self.ps_2D = self.make_2D_ps(height, ps_1D)
        else:
            self.ps_2D = ps_2D
        
        self.dP = self.ps_1D[-1]- self.ps_1D[0]
        
        
    def make_2D_ps(self,height,ps): # sets p(x,y) = p(x) 
        ps_2D = np.zeros((height.Ny, height.Nx))
        
        for i in range(height.Nx):
            for j in range(height.Ny):
                
                y = height.ys[j]
                if y <= height.hs[i]:
                    ps_2D[j,i] = ps[i]
                else:
                    ps_2D[j,i] = None
        # ps_2D = np.flip(ps_2D, axis=0)
        return ps_2D
                    
class FinDiffReynPressure(Pressure):
    def __init__(self, height):
        rhs = fd.make_rhs(height)
        mat = fd.make_mat(height)
        ps_1D = np_fd_solve(mat, rhs)
        super().__init__(height, ps_1D)

class PWLAnalyticReynPressure(Pressure):
    def __init__(self, height):
            
        if not isinstance(height, PWL_Height):
            raise TypeError('Example is not piecewise linear')
        
        rhs = pwl.make_rhs(height)
        linOp = pwl.pwlLinOp(height)
        coefs, exit_code = sp_pwl_gmres(linOp, rhs, tol=1e-10)
         
        if exit_code != 0:
            raise Exception('gmres did not converge')
            
        ps_1D, flux = pwl.make_ps(height, coefs)
        super().__init__(height, ps_1D)
    
class AdjReynPressure(Pressure):
    
    def __init__(self, height):
        reyn_pressure = FinDiffReynPressure(height)
        ps_1D = reyn_pressure.ps_1D
        
        ps_2D = reyn_pressure_adjusted.make_adj_ps(height, ps_1D)

        super().__init__(height, ps_1D, ps_2D)

    