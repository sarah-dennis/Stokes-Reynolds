#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 25 18:24:34 2023

@author: sarahdennis
"""
import pressure_sqrWave as squareWave
import pressure_sawtooth as sawtooth
import pressure_finDiff as finDiff
import pressure_analytic as analytic

class Pressure:
    def __init__(self, p_solver, p_str):
        self.p_solver = p_solver
        self.p_str = p_str
        
    def solve(self):
        return self.p_solver.solve() # = ps
    
class P_Solver:
    def __init__(self, height, p0, pN, solveFun):
        self.height = height
        self.p0 = p0
        self.pN = pN
        self.solve = solveFun
        
#-----------------------------------------------------------------------------
class FinDiff(Pressure):
    def __init__(self, height, p0, pN):
        p_str = "Finite Difference Reynolds ($N_x = %d$)"%height.Nx
        solver = finDiff.solve
        super().__init__(height, p0, pN, solver, p_str)
        
#-----------------------------------------------------------------------------       
class Analytic_Sinusoidal(Pressure):
    def __init__(self, height, p0, pN):
        p_str = "Reynolds Analytic"
        solver = analytic.Solver_Sinusoidal(height, p0, pN)
        super().__init__(height, p0, pN, solver, p_str)

class Analytic_Constant(Pressure):
    def __init__(self, height, p0, pN):
        p_str = " Reynolds Analytic"
        p_solver = analytic.Solver_Constant(height, p0, pN)
        super().__init__(p_solver, p_str)

class Analytic_Linear(Pressure):
    def __init__(self, height, p0, pN):
        p_str = " Reynolds Analytic"
        p_solver = analytic.Solver_Linear(height, p0, pN)
        super().__init__(p_solver, p_str)

class Analytic_Step(Pressure):
    def __init__(self, height, p0, pN):
        p_str = "Reynolds Analytic"
        solver = analytic.stepHeight_solve
        super().__init__(height, p0, pN, solver, p_str)        

#-----------------------------------------------------------------------------        
class PWA_Sawtooth(Pressure):
    def __init__(self, height, p0, pN):
        p_str = "Reynolds Piecewise Analytic"
        solver = sawtooth.Sawtooth_Solver(height, p0, pN)
        super().__init__(height, p0, pN, solver, p_str)
        
#-----------------------------------------------------------------------------
class PWA_RectWave_numpy(Pressure):
    def __init__(self, height, p0, pN):
        p_str = "Reynolds Piecewise Analytic (numpy.linalg)"
        solver = squareWave.numpy_solve(height, p0, pN)
        super().__init__(height, p0, pN, solver, p_str)
        
        
class PWA_RectWave_schurInv(Pressure):
    def __init__(self, height, p0, pN):
        p_str = "Reynolds Piecewise Analytic (schur inv)"
        solver = squareWave.schurInv_solve(height, p0, pN)
        super().__init__(height, p0, pN, solver, p_str)
        
# Current best method ... 
class PWA_RectWave_schurLU(Pressure):
    def __init__(self, height, p0, pN):
        p_str = "Reynolds Piecewise Analytic (schur-LU)"
        solver = squareWave.schurLU_solve(height, p0, pN)
        super().__init__(height, p0, pN, solver, p_str)
