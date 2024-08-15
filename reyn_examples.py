# -*- coding: utf-8 -*-
"""
Created on Wed Feb 22 10:01:42 2023

@author: sarah
"""

import numpy as np

import reyn_heights as heights

import reyn_pressures_pwlinear as p_pwlinear
import reyn_pressures_finDiff as p_finDiff
import reyn_pressures_analytic as p_analytic

import reyn_velocity 
import graphics

class ReynoldsExample:
    def __init__(self, pSolver):
        self.pSolver = pSolver
        self.ps = pSolver.solve()
        self.vel = reyn_velocity.Velocity(self.pSolver.height, self.ps)
        
    def plot_p(self):
        p_title = "Pressure: \n %s"%(self.pSolver.p_str)
        p_labels = ["$x$", "Pressure $p(x)$"]
        graphics.plot_2D(self.ps, self.pSolver.height.xs, p_title, p_labels, color='r')

    def plot_h(self):
        h_title = "Height: \n %s"%(self.pSolver.height.h_str)
        h_labels = ["$x", "Height $h(x)$"]
        graphics.plot_2D(self.pSolver.height.hs, self.pSolver.height.xs, h_title, h_labels, color='r')
    def plot_ph(self):
        ph_title = "Pressure & Height: \n %s"%(self.pSolver.p_str)
        ph_labels = ["Pressure $p(x)$", "Height $h(x)$", "$x$"]
        graphics.plot_2D_twin(self.ps, self.pSolver.height.hs, self.pSolver.height.xs, ph_title, ph_labels)

    def plot_v(self):
        v_title = "Velocity: \n%s"%(self.pSolver.p_str)
        v_ax_labels =  ['$x$', '$y$']
        graphics.plot_stream_height(self.vel.vx, self.vel.vy, self.pSolver.height.hs, self.pSolver.height.xs, self.pSolver.height.ys, v_title, v_ax_labels)
        
#-------------------------------------------------------------------------
# I. Finite Difference
#-------------------------------------------------------------------------
class FinDiff(ReynoldsExample):
    def __init__(self, height, p0, pN):
        pSolver = p_finDiff.Solver_finDiff(height, p0, pN)
        super().__init__(pSolver)

class FinDiff_Rand(FinDiff):
    def __init__(self):
        N = 50
        p0 = 0 
        pN = 10
        x0 = 0
        xf = 1
        h_min = 0.75
        h_max = 1
        U=1
        height = heights.RandomHeight(x0, xf, N, h_min, h_max, U)
        super().__init__(height, p0, pN)

class FinDiff_Custom(FinDiff):
    def __init__(self, example):
        height = example.pSolver.height
        p0 = example.pSolver.p0
        pN = example.pSolver.pN

        super().__init__(height, p0, pN)
# -----------------------------------------------------------------------------
# III. Constant Height
# -----------------------------------------------------------------------------
class Analytic_Constant(ReynoldsExample):
    def __init__(self, height, p0, pN):
        pSolver = p_analytic.Solver_Constant(height, p0, pN)
        super().__init__(pSolver)
        
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
        pSolver = p_analytic.Solver_Linear(height, p0, pN)
        super().__init__(pSolver)
        
class Linear_Ex1(Analytic_Linear):
    def __init__(self):
        N = 100
        U=1
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
        pSolver = p_analytic.Solver_Step(height, p0, pN)
        super().__init__(pSolver)
      
        #TODO testing against Stokes   
class Step_Ex1(Analytic_Step):
    def __init__(self):
        N = 10
        x0 = 0
        xf = 4
        x_step = 1
        h0 = 1
        h1 = 2

        U = 0        
        height = heights.StepHeight(x0, xf, N, h0, h1, x_step, U)
        flux = 1
        H = h1-h0
        p0, pN = reyn_velocity.flux_to_pressure(flux, U,H)
        super().__init__(height, p0, pN)

# -----------------------------------------------------------------------------
# VI. Piecewise Linear
# -----------------------------------------------------------------------------

class PWA_Linear(ReynoldsExample):
    def __init__(self, height, p0, pN):
        pSolver = p_pwlinear.Solver(height, p0, pN) 
        super().__init__(pSolver)

class PiecewiseLinear_Ex1(PWA_Linear):
    def __init__(self):
        N = 500

        x0 = 0
        xf = 1
        
        N_regions = 5
        U=0
        x_peaks = np.array([0, 0.1, 0.5, 0.6, 0.8, 1])

        h_peaks = np.array([[1.0,1.0],[2.0,3.0],[1.0,1.0],[1.0,1.5],[1.0,2.0],[3.0,3.0]])
        height = heights.PiecewiseLinearHeight(x0, xf, N, N_regions, x_peaks, h_peaks,U)
        flux = 1
        H = h_peaks[0]
        p0, pN = reyn_velocity.flux_to_pressure(flux, U, H)
        super().__init__(height, p0, pN)


class PiecewiseLinear_Ex2(PWA_Linear):
    def __init__(self):
        N = 500

        x0 = 0
        xf = 15
        
        N_regions = 10
        U=1
        x_peaks = np.array([0, 1.5, 2.5, 4, 5.2, 7.2, 9, 12, 12.9, 14, 15])

        h_peaks = np.array([[0.5,0.5],[0.3, 0.2],[0.5,0.3],[0.8,0.8],[0.5,0.3],[0.3,0.5],[0.5,0.7],[0.7,0.5],[0.5,0.5],[0.35, 0.25],[0.25,0.25]])
        height = heights.PiecewiseLinearHeight(x0, xf, N, N_regions, x_peaks, h_peaks,U)
        flux = 1
        H = h_peaks[0]
        p0, pN = reyn_velocity.flux_to_pressure(flux, U, H)
        super().__init__(height, p0, pN)
