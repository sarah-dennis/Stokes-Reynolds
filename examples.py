# -*- coding: utf-8 -*-
"""
Created on Wed Feb 22 10:01:42 2023

@author: sarah
"""
import heights as hgt
import pressures as prs
import numpy as np
import domain as dm

#-------------------------------------------------------------------------
# pressure solutions to reynolds equation for the given height functions
#-------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# 0. Constant Height
# -----------------------------------------------------------------------------
def flat(domain, p0, pN, h0):
    height = hgt.ConstantHeight(domain, h0)
    pressure = prs.ConstantPressure(domain, height, p0, pN)
    return height, pressure

# -----------------------------------------------------------------------------
# I. Corrugated Sinusoidal Height
# -----------------------------------------------------------------------------
def corrugated(domain, p0, pN):
    h_mid = 1
    r = 0.5
    k = 2*np.pi
    height = hgt.CorrugatedHeight(domain, h_mid, r, k)
    pressure = prs.CorrugatedPressure(domain, height, p0, pN)
    return height, pressure

# -----------------------------------------------------------------------------
# II. Wedge Height
# -----------------------------------------------------------------------------

def wedge(domain, h0, hf):
    height = hgt.WedgeHeight(domain, h0, hf)
    pressure = prs.WedgePressure(domain, height)
    return height, pressure

def sawtooth(domain, h_min, h_max, n):
    
    h_peaks = np.random.uniform(h_min, h_max, n+1)
    # n = 2
    # h_peaks = [0.001, 2, 0.001]
    Mx = int(domain.Nx//n)
    hs = np.zeros(domain.Nx)
    ps = np.zeros(domain.Nx)
    
    for i in range(n):
        xi_0 = domain.xs[i*Mx]
        xi_f = domain.xs[(i+1)*Mx-1]
        subDomain = dm.Domain(xi_0, xi_f, domain.U, Mx, domain.BC)

        subHeight = hgt.WedgeHeight(subDomain, h_peaks[i], h_peaks[i+1])
        hs[i*Mx : (i+1)*Mx] = subHeight.hs
        
        subPressure = prs.WedgePressure(subDomain, subHeight)
        ps[i*Mx : (i+1)*Mx] = subPressure.ps

    height = hgt.SawtoothHeight(domain, hs)
    pressure = prs.SawtoothPressure(domain, ps)

    return height, pressure

def sawtooth_finDiff(domain, h_min, h_max, n):
    h_peaks = np.random.uniform(h_min, h_max, n+1)
    # h_peaks = [0.001, 2, 0.001]
    Mx = int(domain.Nx//n)
    hs = np.zeros(domain.Nx)

    for i in range(n):
        xi_0 = domain.xs[i*Mx]
        xi_f = domain.xs[(i+1)*Mx-1]
        subDomain = dm.Domain(xi_0, xi_f, domain.U, Mx, domain.BC)
        
        subHeight = hgt.WedgeHeight(subDomain, h_peaks[i], h_peaks[i+1])
        hs[i*Mx : (i+1)*Mx] = subHeight.hs

    height = hgt.SawtoothHeight(domain, hs)
    pressure = prs.FinDiffPressure(domain, height, 0, 0)
    return height, pressure
    

# -----------------------------------------------------------------------------
# III a. Step Height Example (N = 1)
# -----------------------------------------------------------------------------
def step(domain, p0, pN, x_step, h1, h2):

    height = hgt.StepHeight(domain, x_step, h1, h2)
    pressure = prs.StepPressure(domain, height, p0, pN)

    return height, pressure

def variableStepLen_schurLUSolve(domain, p0, pN, x_step, h1, h2):
    height = hgt.StepHeight(domain, x_step, h1, h2)
    if height.n_steps == 1:
        pressure = prs.StepPressure(domain, height, p0, pN)
    else: #N > 1
        pressure = prs.SquareWavePressure_schurLUSolve(domain, height, p0, pN)
    return height, pressure

# -----------------------------------------------------------------------------
# III b. N-Step Height Example (N > 1)
# -----------------------------------------------------------------------------

def squareWave_schurLUSolve(domain, p0, pN, n_steps=25, r=0.001, h_avg=0.1):
    print("\n Loading %d-step Square Wave \n" % (n_steps))
    height = hgt.SquareWaveHeight(domain, h_avg, r, n_steps)
    pressure = prs.SquareWavePressure_schurLUSolve(domain, height, p0, pN)
    return height, pressure


def mySteps_schurLUSolve(domain, p0, pN, h_steps):
    n_steps = len(h_steps)-1
    print("\n Loading %d-step  Wave \n" % (n_steps))
    h_str = " %d-step Wave \n" % (n_steps)
    h_eq =  "h(x) = [...]"
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



