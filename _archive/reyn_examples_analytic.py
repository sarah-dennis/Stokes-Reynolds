# -*- coding: utf-8 -*-
"""
Created on Tue Oct  8 13:43:25 2024

@author: sarah
"""

import reyn_heights as heights

import reyn_pressures_analytic as p_analytic


from reyn_examples import ReynoldsExample

        
# -----------------------------------------------------------------------------
# III. Constant Height
# -----------------------------------------------------------------------------
class Analytic_Constant(ReynoldsExample):
    def __init__(self, height, p0, pN):
        solver = p_analytic.Solver_Constant(height, p0, pN)
        super().__init__(solver)
        
class Constant_Ex1(Analytic_Constant):
    def __init__(self):
        N = 100
        p0 = 0
        pN = 1
        x0 = 0
        xf = 1
        h0 = 1
        height = heights.ConstantHeight(x0, xf, N, h0)
        super().__init__(height, p0, pN)

# -----------------------------------------------------------------------------
# IV. Linear
# -----------------------------------------------------------------------------
class Analytic_Linear(ReynoldsExample):
    def __init__(self, height, p0, pN):
        solver = p_analytic.Solver_Linear(height, p0, pN)
        super().__init__(solver)
        
class Linear_Ex1(Analytic_Linear):
    def __init__(self):
        N = 100
        U = 1
        p0 = 0
        pN = 0
        x0 = 0
        xf = 1
        h0 = 0.2
        h1 = 0.1

        height = heights.LinearHeight(x0, xf, N,h0,h1, U)
        super().__init__(height, p0, pN) 

# -----------------------------------------------------------------------------
# V. Step Height 
# -----------------------------------------------------------------------------
class Analytic_Step(ReynoldsExample):
    def __init__(self, height, p0, pN):
        solver = p_analytic.Solver_Step(height, p0, pN)
        super().__init__(solver)

class Step_Ex1(Analytic_Step):
    def __init__(self):
        N = 100
        x0 = 0
        xf = 4
        x_step = 1
        h0 = 1
        h1 = 2
        U = 0        
        height = heights.StepHeight(x0, xf, N, h0, h1, x_step, U)
        p0 = 1 
        pN = 0
        super().__init__(height, p0, pN)