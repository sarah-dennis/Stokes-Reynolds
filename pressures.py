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
    def __init__(self, height, p0, pN, p_solver, p_str):
        self.height = height
        self.p0 = p0
        self.pN = pN
        self.solve = p_solver # ps = pressure.solve(height, p0, pN)
        self.p_str = p_str

class AnalyticPressure_Linear(Pressure):
    def __init__(self, height, p0, pN):
        p_str = "Analytic Reynolds"
        solver = analytic.linearHeight_solve
        super().__init__(height, p0, pN, solver, p_str)

class AnalyticPressure_Sinusoidal(Pressure): #sinusoidal
    def __init__(self, domain, height, p0, pN):
        p_str = "Analytic Reynolds"
        solver = analytic.sinusoidalHeight_solve
        super().__init__(height, p0, pN, solver, p_str)

class AnalyticPressure_Step(Pressure):
    def __init__(self, domain, height, p0, pN):
        p_str = "Reynolds Analytic"
        solver = analytic.stepHeight_solve
        super().__init__(height, p0, pN, solver, p_str)        

class FinDiffPressure(Pressure):
    def __init__(self, height, p0, pN):
        p_str = "Finite Difference Reynolds ($N_x = %d$)"%height.Nx
        solver = finDiff.solve
        super().__init__(height, p0, pN, solver, p_str)
        
class SawtoothPressure(Pressure):
    def __init__(self, height, p0, pN):
        p_str = "Piecewise Analytic Reynolds "
        solver = sawtooth.solve
        super().__init__(height, p0, pN, solver, p_str)
        

class SquareWavePressure_numpySolve(Pressure):
    def __init__(self, height, p0, pN):
        p_str = "Reynolds Piecewise Analytic (numpy.linalg)"
        solver = squareWave.numpy_solve(height, p0, pN)
        super().__init__(height, p0, pN, solver, p_str)
        
        
class SquareWavePressure_schurInvSolve(Pressure):
    def __init__(self, height, p0, pN):
        p_str = "Reynolds Piecewise Analytic (schur inv)"
        solver = squareWave.schurInv_solve(height, p0, pN)
        super().__init__(height, p0, pN, solver, p_str)
        
# Current best method ... 
class SquareWavePressure_schurLUSolve(Pressure):
    def __init__(self, domain, height, p0, pN):
        p_str = "Reynolds Piecewise Analytic (schur-LU)"
        solver = squareWave.schurLU_solve(height, p0, pN)
        super().__init__(height, p0, pN, solver, p_str)
 