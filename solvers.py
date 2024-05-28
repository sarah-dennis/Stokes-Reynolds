# -*- coding: utf-8 -*-
"""
Created on Tue May 21 16:05:05 2024

@author: sarah
"""

class P_Solver: #TODO rename Reyn_Solver
    def __init__(self, height, p0, pN, solveFun, p_str):
        self.height = height
        self.p0 = p0
        self.pN = pN
        self.solve = solveFun # returns ps #TODO: return velocity??
        self.p_str = p_str
        
class StokesSolver:
    def __init__(self, height, flux, solveFun):
        self.height = height
        self.Q = flux
        self.solve = solveFun #returns psis, us, vs