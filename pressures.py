#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 25 18:24:34 2023

@author: sarahdennis
"""

class P_Solver:
    def __init__(self, height, p0, pN, solveFun, p_str):
        self.height = height
        self.p0 = p0
        self.pN = pN
        self.solve = solveFun # returns ps
        self.p_str = p_str
