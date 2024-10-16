# -*- coding: utf-8 -*-
"""
Created on Wed Feb 22 10:01:42 2023

@author: sarah
"""

import numpy as np

import reyn_heights as heights

import reyn_pressures_pwlinear as p_pwlinear
import reyn_pressures_finDiff as p_finDiff

import reyn_velocity 
import graphics

# NOTE: solver is specifc to example
class ReynoldsExample:
    def __init__(self, solver):
        self.solver = solver
        self.ps, q = solver.solve()
        self.vel = reyn_velocity.Velocity(self.solver.height, self.ps)
        
        dm = self.solver.height
        dP = solver.pN-solver.p0
        
        self.paramstr = '$Re=%.2f$, $dP=%.2f$, $Q=%.3f$, $U=%.1f$'%(dm.Re, dP, q, dm.U)
        
        
        
    def plot_p(self):

        p_title = "Pressure $p(x)$: \n" + self.paramstr
        p_labels = ["$x$", "Pressure $p(x)$"]
        graphics.plot_2D(self.ps, self.solver.height.xs, p_title, p_labels, color='r')

    def plot_h(self):
        h_title = "Height: \n " + self.paramstr
        h_labels = ["$x", "Height $h(x)$"]
        graphics.plot_2D(self.solver.height.hs, self.solver.height.xs, h_title, h_labels, color='r')
    
    def plot_v(self):
        #y-axis reverse
        
        vy=np.flip(self.vel.vy, 0)
        vx=np.flip(self.vel.vx, 0)
        v_title = "Velocity $(u,v)$: \n" + self.paramstr
        v_ax_labels =  ['$|(u,v)|_2$','$x$', '$y$'] 
        uv_mag = np.sqrt(vx**2 + vy**2)
        vmax = 2.5*np.median(uv_mag) #
        graphics.plot_stream_heat(vx, vy, self.solver.height.xs, self.solver.height.ys, uv_mag,v_title, v_ax_labels, vmin=0, vmax=vmax)
        
#-------------------------------------------------------------------------
# I. Discrete (finite difference)
#-------------------------------------------------------------------------

class Discrete_FinDiff(ReynoldsExample): # use height from another example
    def __init__(self, reyn_example):
        height = reyn_example.solver.height
        p0 = reyn_example.solver.p0
        pN = reyn_example.solver.pN
        solver = p_finDiff.Solver_finDiff(height, p0, pN)
        super().__init__(solver)
    
# -----------------------------------------------------------------------------
# II. Piecewise Linear (Piecewise-Analytic)
# -----------------------------------------------------------------------------

class PWL_PWA(ReynoldsExample):
    def __init__(self, pwl_height, p0, pN):
        
        pwl_example = p_pwlinear.Solver(pwl_height, p0, pN) 
        
        super().__init__(pwl_example)


#-----------------------------------------------------------------------------------------------------------------------------------

class BFS(PWL_PWA):
    def __init__(self, U, dP, H, L):
        N = 200

        x0 = 0
        xf = L
        
        N_regions = 2
        x_peaks = np.asarray([x0, 1, xf],float)
        h_peaks = np.asarray([[1,1],[1,H],[H,H]],float)
        height = heights.PiecewiseLinearHeight(x0, xf, N, N_regions, x_peaks, h_peaks,U)

        p0 = -dP
        pN = 0
        super().__init__(height, p0, pN)

class BFS_smooth(PWL_PWA):
    def __init__(self, U, dP, H, L, xL, yL):
        N = 500

        x0 = 0
        xf = L
        x_reattatch=1 +xL
        y_reattatch=H -yL
        N_regions = 3
        x_peaks = np.asarray([x0, 1, x_reattatch, xf],float)
        h_peaks = np.asarray([[1,1],[1,y_reattatch],[y_reattatch,H],[H,H]],float)
        height = heights.PiecewiseLinearHeight(x0, xf, N, N_regions, x_peaks, h_peaks,U)

        p0 = -dP
        pN = 0
        super().__init__(height, p0, pN)

#-----------------------------------------------------------------------------------------------------------------------------------

class HexSlider(PWL_PWA):
    def __init__(self, U, dP):
        N = 200
        
        x0 = 0
        xf = 4

        N_regions = 4
        
        x_peaks = np.asarray([x0, 1, 2, 3, xf],float)
        
        h_peaks = np.asarray(([[1,1],[1,1.5],[2,2],[1.5,1],[1,1]]),float)
        height = heights.PiecewiseLinearHeight(x0, xf, N, N_regions, x_peaks, h_peaks,U)

        p0 = -dP
        pN = 0
        super().__init__(height, p0, pN)















