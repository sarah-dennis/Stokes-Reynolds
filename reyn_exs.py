# -*- coding: utf-8 -*-
"""
Created on Wed Feb 22 10:01:42 2023

@author: sarah
"""

import numpy as np
import heights 

import pressures
#-------------------------------------------------------------------------
# Pressure solutions to reynolds equation for the given height functions
#-------------------------------------------------------------------------

def fd_random():
    Nx = 100 
    ps = np.zeros(Nx)
    
    p0 = 0 
    pN = 0 
    
    x0 = 0
    xf = 1

    h_min = 0.1 
    h_max = 0.5
    height = heights.Random(x0, xf, h_min, h_max, Nx)

    pressure = pressures.FinDiffPressure(height, p0, pN)
    
    ps = pressure.solve()
    
    return ps

def fd_blank():
    Nx = 100
    ps = np.zeros(Nx)
    
    # p0 = 0
    # pN = 0
    
    # x0 = 0
    # xf = 1

    # height = ...
    # pressure = ...
    # ps = pressure.solve()
    
    return ps

# -----------------------------------------------------------------------------
# I. Sinusoidal Height
# -----------------------------------------------------------------------------
def analytic_sinusoidal():
    Nx = 100
    ps = np.zeros(Nx)
    
    # p0 = 0
    # pN = 0
    
    # x0 = 0
    # xf = 1
    
    # height = ...
    # pressure = ...
    # ps = pressure.solve()

    return ps

# -----------------------------------------------------------------------------
# II. Linear Height
# -----------------------------------------------------------------------------


# II. a) Piecewise-linear -> piecewise solution (from double int. Reynolds )

def sawtooth(domain, p0, pN, h_min, h_max, N): 
    N = 2
    
    #uniform dx
    xi = domain.x0 + np.arange(0, N+1) * (domain.xf - domain.x0)/N

    #I. random height
    # hi = np.random.uniform(h_min, h_max, N+1)
    
    #II. prescribed height
    # hi = np.array([h_min + .3*h_max, h_min + .8*h_max, h_min +.2*h_max, h_min +  .5*h_max, h_min + .1*h_max, h_min + h_max, h_min + .4*h_max])
    # hi = np.array([h_min + .3*h_max, h_min + .4*h_max, h_min +.2*h_max, h_min +  .35*h_max, h_min + .1*h_max, h_min + .4*h_max, h_min + .25*h_max])
    hi = np.array([h_min, h_max, h_min])
    height = hgt.SawtoothHeight(domain, xi, hi)
    
    st_pressure = prs.SawtoothPressure(domain, height, p0, pN)
    


    return height, st_pressure


def sawtoothRand(domain, p0, pN, h_min, h_max, N): 
    
    #random dx
    xi = np.sort(np.random.uniform(domain.x0, domain.xf, N+1))
    xi[0] = domain.x0
    xi[-1] = domain.xf
    
    hi = np.random.uniform(h_min, h_max, N+1)
    
    height = hgt.SawtoothHeight(domain, xi, hi)
    
    pressure = prs.SawtoothPressure(domain, height, p0, pN)



    return height, pressure


# -----------------------------------------------------------------------------
# III a. Step Height Example (N = 1)
# -----------------------------------------------------------------------------
def step(domain, p0, pN, x_step, h1, h2):

    height = hgt.RayleighStepHeight(domain, x_step, h1, h2)
    pressure = prs.RayleighStepPressure(domain, height, p0, pN)

    return height, pressure


# -----------------------------------------------------------------------------
# III b. N-Step Height Example (N > 1)
# -----------------------------------------------------------------------------

def squareWave_schurLUSolve(domain, p0, pN, n_steps=25, r=0.001, h_avg=0.1):
    print("\n Loading %d-step Square Wave \n" % (n_steps))
    height = hgt.SquareWaveHeight(domain, h_avg, r, n_steps)
    pressure = prs.SquareWavePressure_schurLUSolve(domain, height, p0, pN)
    return height, pressure


def randSteps_schurLUSolve(domain, p0, pN, h_min, h_max, n_steps):
    print("\n Loading %d-step  Wave \n" % (n_steps))
    h_str = " %d-step Wave \n" % (n_steps)
    h_eq =  "h(x) \in [%.2f, ..., %.2f]"%(h_min, h_max)
    h_steps = np.random.uniform(h_min, h_max, n_steps+1)
    height = hgt.NStepHeight(domain, n_steps, h_steps, h_str, h_eq)
    pressure = prs.SquareWavePressure_schurLUSolve(domain, height, p0, pN)
    return height, pressure

def squareWave_schurInvSolve(domain, p0, pN, n_steps=205, r=0.001, h_avg=0.1):
    print("\n Loading %d-step Square Wave \n" % (n_steps))
    height = hgt.SquareWaveHeight(domain, h_avg, r, n_steps)
    pressure = prs.SquareWavePressure_schurInvSolve(domain, height, p0, pN)
    return height, pressure


def squareWave_pySolve(domain, p0, pN, n_steps=25, r=0.001, h_avg=0.1):
    print("\n Loading %d-step Square Wave \n" % (n_steps))
    height = hgt.SquareWaveHeight(domain, h_avg, r, n_steps)
    pressure = prs.SquareWavePressure_pySolve(domain, height, p0, pN)
    return height, pressure

def squareWave_schurGmresSolve(domain, p0, pN, n_steps=2105, r=0.001, h_avg=0.1):
    print("\n Loading %d-step Square Wave \n" % (n_steps))
    height = hgt.SquareWaveHeight(domain, h_avg, r, n_steps)
    pressure = prs.SquareWavePressure_schurGmresSolve(domain, height, p0, pN)
    return height, pressure


