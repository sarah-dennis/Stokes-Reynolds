 #!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec  5 13:28:42 2024

@author: sarahdennis
"""
from pressure import Reyn_Pressure
import reyn_pressure_finDiff as fd
import reyn_pressure_pwl as pwl
import reyn_pressure_pwc as pwc

                    
class FinDiff_ReynPressure(Reyn_Pressure):
    def __init__(self, height, BC):
        ps_1D = fd.fd_solve(height, BC)
        
        super().__init__(height, ps_1D)


class PwlGMRes_ReynPressure(Reyn_Pressure):
    def __init__(self, height, BC):
            

            
        ps_1D = pwl.gmres_solve(height, BC)
        super().__init__(height,ps_1D)
        
class PwlSchur_ReynPressure(Reyn_Pressure):
    def __init__(self, height, BC):
            
        ps_1D = pwl.schur_solve(height, BC)
        super().__init__(height,  ps_1D)
    
class PwcSchur_ReynPressure(Reyn_Pressure):
    def __init__(self, height, BC): 
        
        ps_1D = pwc.schur_solve(height, BC)

        super().__init__(height,ps_1D)
        
class PwcSchur_parallel_ReynPressure(Reyn_Pressure):
    def __init__(self, height, BC): 
        
        ps_1D = pwc.schur_solve_parallel(height, BC)


        super().__init__(height, ps_1D)
