# -*- coding: utf-8 -*-
"""
Created on Wed Feb 22 10:01:42 2023

@author: sarah
"""
import heights as hgt
import pressures as prs
import numpy as np

# Analytic pressure solutions to given height functions
#------------------------------------------------------------------------------
# 0. Constant Height 
#------------------------------------------------------------------------------
def flat(domain, p0, pN):
    h0 = 0.5
    #height = hgt.constHeight(domain, h0)
    height = hgt.SquareWaveHeight(domain, h0, r=0, n_steps=1)
    #pressure = exp.UnknownPressure(domain, height, p0, pN)
    pressure = prs.SquareWavePressure(domain, height, p0, pN)
    return height, pressure

#------------------------------------------------------------------------------
# I. Corrugated Sinusoidal Height 
#------------------------------------------------------------------------------
def corrugated(domain, p0, pN):
    h_mid = 1 
    r = 0.5
    k = 2*np.pi
    height = hgt.CorrugatedHeight(domain, h_mid, r, k)
    pressure = prs.CorrugatedPressure(domain, height, p0, pN)
    return height, pressure

#------------------------------------------------------------------------------
# II. Wedge Height 
#------------------------------------------------------------------------------
def wedge(domain, p0, pN):
    h_min = 0.1
    m = -2
    height = hgt.WedgeHeight(domain, h_min, m)
    pressure = prs.WedgePressure(domain, height, p0, pN)

    return height, pressure

#------------------------------------------------------------------------------
# III a. Step Height Example
#------------------------------------------------------------------------------
def step(domain, p0, pN):
    h_left = 0.3
    h_right = 0.1
    height = hgt.StepHeight(domain, h_left, h_right)
    pressure = prs.StepPressure(domain, height, p0, pN)
    return height, pressure

#------------------------------------------------------------------------------
# III b. Two Step Height Example
#------------------------------------------------------------------------------
def twoStep(domain, p0, pN):
    h_left = 2
    h_center = 1
    h_right = 2
    height = hgt.TwoStepHeight(domain, h_left, h_center, h_right)
    pressure = prs.TwoStepPressure(domain, height, p0, pN)
    return height, pressure

#------------------------------------------------------------------------------
# III c. N-Step Height Example
#------------------------------------------------------------------------------
def squareWave_schurLUSolve(domain, p0, pN, n_steps=205, r=0.001, h_avg=0.1):
    if domain.Nx < n_steps * 3:
        print("Warning: Nx < nsteps * 3")   
    print("Loading %d-step Square Wave \n"%(n_steps))
    height = hgt.SquareWaveHeight(domain, h_avg, r, n_steps)
    pressure = prs.SquareWavePressure_schurLUSolve(domain, height, p0, pN)

    return height, pressure

def squareWave_schurInvSolve(domain, p0, pN, n_steps=205, r=0.001, h_avg=0.1):
    if domain.Nx < n_steps * 3:
        print("Warning: Nx < nsteps * 3")   
    print("Loading %d-step Square Wave \n"%(n_steps))
    height = hgt.SquareWaveHeight(domain, h_avg, r, n_steps)
    pressure = prs.SquareWavePressure_schurInvSolve(domain, height, p0, pN)

    return height, pressure

def squareWave_pySolve(domain, p0, pN, n_steps=205, r=0.001, h_avg=1):
    if domain.Nx < n_steps * 3:
        print("Warning: Nx < nsteps * 3")
        
    print("Loading %d-step Square Wave \n"%(n_steps))
    height = hgt.SquareWaveHeight(domain, h_avg, r, n_steps)
    pressure = prs.SquareWavePressure_pySolve(domain, height, p0, pN)

    return height, pressure