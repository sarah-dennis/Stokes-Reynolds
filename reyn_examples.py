# -*- coding: utf-8 -*-
"""
Created on Wed Feb 22 10:01:42 2023

@author: sarah
"""

import numpy as np

import reyn_heights as heights

import reyn_pressures_stepWave as p_stepWave
import reyn_pressures_sawtooth as p_sawtooth
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
        # graphics.plot_quiver_height(self.vel.vx, self.vel.vy, self.pSolver.height.hs, self.pSolver.height.xs, self.pSolver.height.ys, v_title, v_ax_labels)
    
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
        height = heights.RandomHeight(x0, xf, N, h_min, h_max)
        super().__init__(height, p0, pN)

class FinDiff_Custom(FinDiff):
    def __init__(self, example):
        height = example.height
        p0 = example.p0
        pN = example.pN

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
        pSolver = p_sawtooth.Solver(height, p0, pN)
        super().__init__(pSolver)
        
class Linear_Ex1(Analytic_Linear):
    def __init__(self):
        N = 100
        p0 = 0
        pN = 100
        x0 = 0
        xf = 1
        h0 = 0.1
        h1 = 0.1
        height = heights.SawtoothHeight(x0, xf, N, 1, [x0,xf], [h0,h1])
        super().__init__(height, p0, pN) 


# -----------------------------------------------------------------------------
# V. Step Height 
# -----------------------------------------------------------------------------
class Analytic_Step(ReynoldsExample):
    def __init__(self, height, p0, pN):
        pSolver = p_analytic.Solver_Step(height, p0, pN)
        self.height = height
        self.p0 = p0
        self.pN = pN
        super().__init__(pSolver)
        
class Step_Ex1(Analytic_Step):
    def __init__(self):
        N = 100
        p0 = 0
        pN = 0
        x0 = 0
        xf = 4
        x_step = (xf - x0) / 4
        h0 = 1
        h1 = 2
        height = heights.StepHeight(x0, xf, N, h0, h1, x_step)
        super().__init__(height, p0, pN)

# -----------------------------------------------------------------------------
# VI. Step Wave **
# -----------------------------------------------------------------------------
class PWA_StepWave(ReynoldsExample):
    def __init__(self, height, p0, pN):
        pSolver = p_stepWave.Solver_schurLU(height, p0, pN)
        # pSolver = p_stepWave.Solver_schurInv(height, p0, pN)
        # pSolver = p_stepWave.Solver_numpy(height, p0, pN)
        self.height=height
        super().__init__(pSolver)
        
class StepWave_Ex1(PWA_StepWave):
    def __init__(self):
        N = 1000 #this is just used for plotting
        p0 = 0
        pN = 0
        x0 = 0
        xf = 1
        N_steps = 5
        h_min = 0.1 
        h_max = 0.5
        h_steps = np.zeros(N_steps+1)

        # symmetric oscillation
        h_avg = (h_max + h_min)/2
        r = (h_max - h_min)/2
        for i in range(N_steps+1):
                h_steps[i] = h_avg + (-1)**i * r
        
        height = heights.StepWaveHeight(x0, xf, N, N_steps, h_steps)
        super().__init__(height, p0, pN)

class StepWave_Ex2(PWA_StepWave):
    def __init__(self):
        N = 1000 #this is just used for plotting
        p0 = 0
        pN = 0
        x0 = 0
        xf = 1
        N_steps = 3
        h_min = .1
        h_max = 0.2
        h_steps = np.zeros(N_steps+1)
        
        # random wave
        h_steps = np.random.uniform(h_min, h_max, N_steps+1)
        
        height = heights.StepWaveHeight(x0, xf, N, N_steps, h_steps)
        super().__init__(height, p0, pN)


class StepWave_Ex3(PWA_StepWave):
    def __init__(self):
        N = 200 #this is just used for plotting
        p0 = 0
        pN = 0
        x0 = 0
        xf = 4
        N_steps = 3
        h_min = 1/(3*N)
        
        h_steps = np.array([h_min, 0.5, 1, h_min])
        
        # random wave

        height = heights.StepWaveHeight(x0, xf, N, N_steps, h_steps)
        super().__init__(height, p0, pN)

# -----------------------------------------------------------------------------
# VI. Sawtooth **
# -----------------------------------------------------------------------------

class PWA_Sawtooth(ReynoldsExample):
    def __init__(self, height, p0, pN):
        pSolver = p_sawtooth.Solver(height, p0, pN)
        super().__init__(pSolver)
        
class Sawtooth_Ex1(PWA_Sawtooth):
    def __init__(self):
        N = 100
        p0 = 0
        pN = 0
        x0 = 0
        xf = 1
        
        N_regions = 5

        x_peaks = np.array([0, 0.1, 0.5, 0.6, 0.8, 1])

        h_peaks = np.array([1, 2, 1, 2, 1, 2])
        height = heights.SawtoothHeight(x0, xf, N, N_regions, x_peaks, h_peaks)
        super().__init__(height, p0, pN)

class Sawtooth_Ex2(PWA_Sawtooth):
    def __init__(self):
        N = 100
        p0 = 0
        pN = 0
        x0 = 0
        xf = 10
        
        N_regions = 8
        h_min = 0.1
        h_max =0.3
        
        # uniform width 
        x_peaks = x0 + np.arange(0, N_regions+1) * (xf - x0)/N_regions

        # random height
        h_peaks = np.random.uniform(h_min, h_max, N_regions+1)

        height = heights.SawtoothHeight(x0, xf, N, N_regions, x_peaks, h_peaks)
        super().__init__(height, p0, pN)
        
class Sawtooth_Ex3(PWA_Sawtooth):
    def __init__(self):
        N = 100
        p0 = 0
        pN = 0
        x0 = 0
        xf = 10
        
        N_regions = 8
        h_min = 0.1
        h_max =0.3
        
        # uniform width 
        x_peaks = x0 + np.arange(0, N_regions+1) * (xf - x0)/N_regions

        # random height
        h_peaks = np.random.uniform(h_min, h_max, N_regions+1)
        h_peaks[4] = h_peaks[3]
        height = heights.SawtoothHeight(x0, xf, N, N_regions, x_peaks, h_peaks)
        super().__init__(height, p0, pN)
        
     
# -----------------------------------------------------------------------------
# VI. Sawtooth **
# -----------------------------------------------------------------------------

class PWA_Linear(ReynoldsExample):
    def __init__(self, height, p0, pN):
        pSolver = p_pwlinear.Solver(height, p0, pN) 
        self.height=height
        self.p0 = p0
        self.pN = pN
        super().__init__(pSolver)

class PiecewiseLinear_Ex0(PWA_Linear):
    def __init__(self):
        N = 500
        p0 = 0
        pN = 0
        x0 = 0
        xf = 1
        
        N_regions = 5

        x_peaks = np.array([0, 0.1, 0.5, 0.6, 0.65, 1])

        h_peaks = np.array([[1,1], [2,3], [1,1], [2,2], [1,3], [1,1]])
        height = heights.PiecewiseLinearHeight(x0, xf, N, N_regions, x_peaks, h_peaks)
        super().__init__(height, p0, pN)

class PiecewiseLinear_Ex1(PWA_Linear):
    def __init__(self):
        N = 500
        p0 = 0
        pN = 0
        x0 = 0
        xf = 1
        
        N_regions = 5

        x_peaks = np.array([0, 0.1, 0.5, 0.6, 0.8, 1])

        h_peaks = np.array([[1.0,1.0],[2.0,3.0],[1.0,1.0],[1.0,1.5],[1.0,2.0],[3.0,3.0]])
        height = heights.PiecewiseLinearHeight(x0, xf, N, N_regions, x_peaks, h_peaks)
        super().__init__(height, p0, pN)


class PiecewiseLinear_Ex2(PWA_Linear):
    def __init__(self):
        N = 500
        p0 = 1
        pN = 0
        x0 = 0
        xf = 15
        
        N_regions = 10

        x_peaks = np.array([0, 1.5, 2.5, 4, 5.2, 7.2, 9, 12, 12.9, 14, 15])

        h_peaks = np.array([[0.5,0.5],[0.3, 0.2],[0.5,0.3],[0.8,0.8],[0.5,0.3],[0.3,0.5],[0.5,0.7],[0.7,0.5],[0.5,0.5],[0.35, 0.25],[0.25,0.25]])
        height = heights.PiecewiseLinearHeight(x0, xf, N, N_regions, x_peaks, h_peaks)
        super().__init__(height, p0, pN)
