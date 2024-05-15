# -*- coding: utf-8 -*-
"""
Created on Wed Feb 22 10:01:42 2023

@author: sarah
"""

import numpy as np

import heights 

import pressures_stepWave as p_stepWave
import pressures_sawtooth as p_sawtooth
import pressures_finDiff as p_finDiff
import pressures_analytic as p_analytic

import velocity 
import graphics

class ReynoldsExample:
    def __init__(self, pSolver):
        self.pSolver = pSolver
        self.ps = pSolver.solve()
        self.vel = velocity.Velocity(self.pSolver.height, self.ps)
        
    def plot_p(self):
        p_title = "Pressure: \n %s"%(self.pSolver.p_str)
        p_labels = ["$x$", "Pressure $p(x)$", ]
        graphics.plot_2D(self.ps, self.pSolver.height.xs, p_title, p_labels)

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

class FinDiff_Ex1(FinDiff):
    def __init__(self):
        Nx = 100 
        p0 = 0 
        pN = 0 
        x0 = 0
        xf = 1
        h_min = 0.1 
        h_max = 0.15
        height = heights.Random(x0, xf, h_min, h_max, Nx)
        super().__init__(height, p0, pN)

# -----------------------------------------------------------------------------
# II. Sinusoidal Height
# -----------------------------------------------------------------------------
class Analytic_Sinusoidal(ReynoldsExample):
    def __init__(self, height, p0, pN):
        pSolver = p_analytic.Solver_Sinusoidal(height, p0, pN)
        super().__init__(pSolver)
        
class Sinusoidal_Ex1(Analytic_Sinusoidal):
    def __init__(self):
        Nx = 100
        p0 = 0
        pN = 0
        x0 = 0
        xf = 10
        h_avg = 0.5
        r = 0.2 
        k = 2
        height = heights.SinsusoidalHeight(x0, xf, Nx, h_avg, r, k)
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
        Nx = 100
        p0 = 0
        pN = 0
        x0 = 0
        xf = 1
        h0 = 1
        height = heights.ConstantHeight(x0, xf, Nx, h0)
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
        Nx = 100
        p0 = 0
        pN = 0
        x0 = 0
        xf = 1
        h0 = 0.1
        h1 = 0.01 # h1 != h0
        height = heights.LinearHeight(x0, xf, Nx, h0, h1)
        super().__init__(height, p0, pN) 


# -----------------------------------------------------------------------------
# V. Step Height 
# -----------------------------------------------------------------------------
class Analytic_Step(ReynoldsExample):
    def __init__(self, height, p0, pN):
        pSolver = p_analytic.Solver_Step(height, p0, pN)
        super().__init__(pSolver)
        
class Step_Ex1(Analytic_Step):
    def __init__(self):
        Nx = 100
        p0 = 0
        pN = 0
        x0 = 0
        xf = 1
        x_step = (xf - x0) / 3
        h0 = 0.1
        h1 = 0.1
        height = heights.StepHeight(x0, xf, Nx, h0, h1, x_step)
        super().__init__(height, p0, pN)

# -----------------------------------------------------------------------------
# VI. Step Wave **
# -----------------------------------------------------------------------------
class PWA_StepWave(ReynoldsExample):
    def __init__(self, height, p0, pN):
        # pSolver = p_stepWave.Solver_schurLU(height, p0, pN)
        # pSolver = p_stepWave.Solver_schurInv(height, p0, pN)
        pSolver = p_stepWave.Solver_numpy(height, p0, pN)
        super().__init__(pSolver)
        
class StepWave_Ex1(PWA_StepWave):
    def __init__(self):
        Nx = 100
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
        
        height = heights.StepWaveHeight(x0, xf, Nx, N_steps, h_steps)
        super().__init__(height, p0, pN)

class StepWave_Ex2(PWA_StepWave):
    def __init__(self):
        Nx = 100
        p0 = 0
        pN = 0
        x0 = 0
        xf = 1
        N_steps = 5
        h_min = 0.1 
        h_max = 0.5
        h_steps = np.zeros(N_steps+1)
        
        # random wave
        h_steps = np.random.uniform(h_min, h_max, N_steps+1)
        
        height = heights.StepWaveHeight(x0, xf, Nx, N_steps, h_steps)
        super().__init__(height, p0, pN)

# -----------------------------------------------------------------------------
# VI. Sawtooth **
# -----------------------------------------------------------------------------
#TODO: does PWA_Sawtooth break for an interval of constant height?
class PWA_Sawtooth(ReynoldsExample):
    def __init__(self, height, p0, pN):
        pSolver = p_sawtooth.Solver(height, p0, pN)
        super().__init__(pSolver)
        
class Sawtooth_Ex1(PWA_Sawtooth):
    def __init__(self):
        Nx = 100
        p0 = 0
        pN = 0
        x0 = 0
        xf = 1
        
        N_regions = 5
        h_min = 0.1 
        h_max = 0.5
        
        # uniform width 
        # x_peaks = x0 + np.arange(0, N_regions+1) * (xf - x0)/N_regions
        x_peaks = np.array([0, 0.1, 0.5, 0.6, 0.8, 1])

        # random height
        # h_peaks = np.random.uniform(h_min, h_max, N_regions+1)
        h_peaks = np.array([1, 2, 1, 2, 1, 2])
        height = heights.SawtoothHeight(x0, xf, Nx, N_regions, x_peaks, h_peaks)
        super().__init__(height, p0, pN)

