# -*- coding: utf-8 -*-
"""
Created on Tue Oct  1 11:56:12 2024

link to bash file: sd_run.sh
@author: sarah
"""

import stokes_control 
import stokes_examples as examples

stokes_control.example = examples.HexSlider_Re05_Q2
stokes_control.new_run_many(10,2,6)
