# -*- coding: utf-8 -*-
"""
Created on Tue May 21 16:05:05 2024

@author: sarah
"""

class Reynolds_Solver: 
    def __init__(self, height, p0, pN, solveFun, p_str):
        self.height = height
        self.p0 = p0
        self.pN = pN
        self.solve = solveFun # solveFun(): height -> ps
        self.filestr = p_str
