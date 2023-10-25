# -*- coding: utf-8 -*-
"""
Created on Wed Feb 22 10:01:42 2023

@author: sarah
"""
import heights as hgt
import pressures as prs
import velocity as vel
import numpy as np


# pressure solutions to reynolds equation for the given height functions

# ------------------------------------------------------------------------------
# 0. Constant Height
# ------------------------------------------------------------------------------
def flat(domain, p0, pN, h0):
    height = hgt.ConstantHeight(domain, h0)
    pressure = prs.ConstantPressure(domain, height, p0, pN)
    return height, pressure

# ------------------------------------------------------------------------------
# I. Corrugated Sinusoidal Height
# ------------------------------------------------------------------------------


def corrugated(domain, p0, pN):
    h_mid = 1
    r = 0.5
    k = 2*np.pi
    height = hgt.CorrugatedHeight(domain, h_mid, r, k)
    pressure = prs.CorrugatedPressure(domain, height, p0, pN)
    return height, pressure

# ------------------------------------------------------------------------------
# II. Wedge Height
# ------------------------------------------------------------------------------


def wedge(domain, p0, pN):
    h_min = 0.1
    m = -2
    height = hgt.WedgeHeight(domain, h_min, m)
    pressure = prs.WedgePressure(domain, height, p0, pN)
    return height, pressure


# ------------------------------------------------------------------------------
# III a. Step Height Example
# ------------------------------------------------------------------------------
def step(domain, p0, pN, x_step, r, h_avg):

    height = hgt.StepHeight(domain, x_step, h_avg, r)
    pressure = prs.StepPressure(domain, height, p0, pN)

    return height, pressure

# ------------------------------------------------------------------------------
# III b. N-Step Height Example
# ------------------------------------------------------------------------------


def squareWave_schurLUSolve(domain, p0, pN, n_steps=25, r=0.001, h_avg=0.1):
    print("\n Loading %d-step Square Wave \n" % (n_steps))
    height = hgt.SquareWaveHeight(domain, h_avg, r, n_steps)
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

def squareWave_gmresSolve(domain, p0, pN, n_steps=2105, r=0.001, h_avg=0.1):
    print("\n Loading %d-step Square Wave \n" % (n_steps))
    height = hgt.SquareWaveHeight(domain, h_avg, r, n_steps)
    pressure = prs.SquareWavePressure_gmresSolve(domain, height, p0, pN)
    return height, pressure


def mySteps_schurLUSolve(domain, p0, pN, h_steps):
    n_steps = len(h_steps)-1
    print("\n Loading %d-step  Wave \n" % (n_steps))
    h_str = " %d-step Wave \n" % (n_steps)
    h_eq =  "h(x) = [...]"
    height = hgt.NStepHeight(domain, n_steps, h_steps, h_str, h_eq)
    pressure = prs.SquareWavePressure_schurLUSolve(domain, height, p0, pN)
    return height, pressure
