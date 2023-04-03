# -*- coding: utf-8 -*-
"""
Created on Wed Feb 22 10:01:42 2023

@author: sarah
"""
import heights as hgt
import exactPressures as pres
import numpy as np

#------------------------------------------------------------------------------
# I. Corrugated Height example
#------------------------------------------------------------------------------
def corrugated(domain, p0, pN):
    h_mid = 1 
    r = 0.5
    k = 2*np.pi
    height = hgt.CorrugatedHeight(domain, h_mid, r, k)
    pressure = pres.CorrugatedPressure(domain, height, p0, pN)
    return height, pressure

#------------------------------------------------------------------------------
# II. Wedge Height example
#------------------------------------------------------------------------------
def wedge(domain, p0, pN):
    h_min = 0.1
    m = -2
    height = hgt.WedgeHeight(domain, h_min, m)
    pressure = pres.WedgePressure(domain, height, p0, pN)

    return height, pressure

#------------------------------------------------------------------------------
# III a. Step Height Example
#------------------------------------------------------------------------------
def step(domain, p0, pN):
    h_left = 0.3
    h_right = 0.1
    height = hgt.StepHeight(domain, h_left, h_right)
    pressure = pres.StepPressure(domain, height, p0, pN)
    return height, pressure

#------------------------------------------------------------------------------
# III b. Two Step Height Example
#------------------------------------------------------------------------------
def twoStep(domain, p0, pN):
    h_left = 2
    h_center = 1
    h_right = 2
    height = hgt.TwoStepHeight(domain, h_left, h_center, h_right)
    pressure = pres.TwoStepPressure(domain, height, p0, pN)
    return height, pressure

#------------------------------------------------------------------------------
# IV c. N-Step Height Example
#------------------------------------------------------------------------------

def squareWave(domain, p0, pN, n_steps=35, r=0.1, h_avg=0.2):
    if domain.Nx < n_steps * 3:
        print("Warning: Nx < nsteps * 3")
    
    height = hgt.SquareWaveHeight(domain, h_avg, r, n_steps)
    pressure = pres.SquareWavePressure(domain, height, p0, pN)
    return height, pressure

#------------------------------------------------------------------------------
# VI e. Const Height Example
#------------------------------------------------------------------------------
def flat(domain, p0, pN):
    h0 = 0.5
    #height = hgt.constHeight(domain, h0)
    height = hgt.SquareWaveHeight(domain, h0, r=0, n_steps=1)
    #pressure = exp.UnknownPressure(domain, height, p0, pN)
    pressure = pres.SquareWavePressure(domain, height, p0, pN)
    return height, pressure

#------------------------------------------------------------------------------
# V d. Dimple Height Example
#------------------------------------------------------------------------------
def dimple(domain, p0, pN):
    h_min = 0.5
    h_max = 1
    len_ratio = 64
    height = hgt.DimpleHeight(domain, len_ratio, h_max, h_min)
    pressure = pres.UnknownPressure(domain, height, p0, pN)
    return height, pressure