# -*- coding: utf-8 -*-
"""
Created on Tue Aug  6 09:17:16 2024

@author: sarah
"""


# -----------------------------------------------------------------------------
# VI. Step Wave  -- redundant by PWL
# -----------------------------------------------------------------------------
class PWA_StepWave(ReynoldsExample):
    def __init__(self, height, p0, pN):
        pSolver = p_stepWave.Solver_schurLU(height, p0, pN)
        # pSolver = p_stepWave.Solver_schurInv(height, p0, pN)
        # pSolver = p_stepWave.Solver_numpy(height, p0, pN)
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
        U=1
        # symmetric oscillation
        h_avg = (h_max + h_min)/2
        r = (h_max - h_min)/2
        for i in range(N_steps+1):
                h_steps[i] = h_avg + (-1)**i * r
        
        height = heights.StepWaveHeight(x0, xf, N, N_steps, h_steps, U)
        super().__init__(height, p0, pN)

class StepWave_Ex2(PWA_StepWave):
    def __init__(self):
        N = 1000 #this is just used for plotting
        p0 = 0
        pN = 0
        x0 = 0
        xf = 1
        U=1
        N_steps = 3
        h_min = .1
        h_max = 0.2
        h_steps = np.zeros(N_steps+1)
        
        # random wave
        h_steps = np.random.uniform(h_min, h_max, N_steps+1)
        
        height = heights.StepWaveHeight(x0, xf, N, N_steps, h_steps, U)
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
        U = 1
        h_steps = np.array([h_min, 0.5, 1, h_min])
        
        # random wave

        height = heights.StepWaveHeight(x0, xf, N, N_steps, h_steps, U)
        super().__init__(height, p0, pN)

# -----------------------------------------------------------------------------
# VI. Sawtooth -- redundant by PWL
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
        
     