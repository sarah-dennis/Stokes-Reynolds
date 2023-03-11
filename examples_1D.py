# -*- coding: utf-8 -*-
"""
Created on Wed Feb 22 10:01:42 2023

@author: sarah
"""
import heights as hgt
import exactPressures as exp
import numpy as np

#------------------------------------------------------------------------------
# I. Corrugated Height example
#------------------------------------------------------------------------------
# Jan 29: converges 2nd order 
def corrugated(domain):
    h_mid = 1 
    r = 0.5
    k = 2*np.pi
    height = hgt.CorrugatedHeight(domain, h_mid, r, k)
    p0 = 0
    pN = 0
    pressure = exp.CorrugatedPressure(domain, height, p0, pN)
    return height, pressure

#------------------------------------------------------------------------------
# II. Wedge Height example
#------------------------------------------------------------------------------
def wedge(domain):
    h_min = 0.1
    m = -2
    height = hgt.WedgeHeight(domain, h_min, m)
    p0 = 0
    pN = 0
    pressure = exp.WedgePressure(domain, height, p0, pN)

    return height, pressure

#------------------------------------------------------------------------------
# III a. Step Height Example
#------------------------------------------------------------------------------
# Jan 29: converges 1st order 
def step(domain):
    h_left = 0.3
    h_right = 0.1
    height = hgt.StepHeight(domain, h_left, h_right)
    p0 = 0
    pN = 0
    pressure = exp.StepPressure(domain, height, p0, pN)
    return height, pressure

#------------------------------------------------------------------------------
# III b. Two Step Height Example
#------------------------------------------------------------------------------
# Jan 29: converges 1st order 
def twoStep(domain):
    h_left = 2
    h_center = 1
    h_right = 2
    height = hgt.TwoStepHeight(domain, h_left, h_center, h_right)
    p0 = 0
    pN = 0
    pressure = exp.TwoStepPressure(domain, height, p0, pN)
    return height, pressure

#------------------------------------------------------------------------------
# IV c. N-Step Height Example
#------------------------------------------------------------------------------
# Mar 6: converges 1st order 
def squareWave(domain, n_steps=25, r=0.1):
    h_avg = 1
    
    if domain.Nx < n_steps * 3:
        print("Warning: Nx < nsteps * 3")
    
    height = hgt.SquareWaveHeight(domain, h_avg, r, n_steps)
    p0 = 0
    pN = 0
    pressure = exp.SquareWavePressure(domain, height, p0, pN)
    return height, pressure

