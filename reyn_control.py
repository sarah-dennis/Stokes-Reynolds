#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 22 15:05:40 2023
@author: sarahdennis
"""
import csv
import reyn_examples as rex

# solution = rex.FinDiff_Ex1()
# solution = rex.FinDiff_Ex2()
# solution = rex.Constant_Ex1()
# solution = rex.Linear_Ex1()       
# solution = rex.Step_Ex1()
# solution = rex.StepWave_Ex1()
# solution = rex.StepWave_Ex2()
# solution = rex.Sawtooth_Ex1()
# solution = rex.Sawtooth_Ex2()
solution = rex.Sawtooth_Ex3()

# solution = rex.StepWave_Ex3()

#------------------------------------------------------------------------------
# plotting 
#------------------------------------------------------------------------------
# solution.plot_ph()
solution.plot_p()
# solution.plot_h()
solution.plot_v()

#------------------------------------------------------------------------------
# Error
#------------------------------------------------------------------------------

# infNorm_err = np.max(np.abs(al_pressure.ps - fd_pressure.ps))
# print("Analytic to Numerical Error: %.8f"%infNorm_err)

# num_err_title = "Numerical Error for %s"%height.h_str
# num_err_axis = ["$x$", "Error (Analytic - Numerical)"]
# graph.plot_2D(al_pressure.ps - fd_pressure.ps, domain.xs, num_err_title, num_err_axis)

#------------------------------------------------------------------------------
# CSV writing
#------------------------------------------------------------------------------

def write_solution(filename, length, fd_ps, al_ps, N, err):

    with open(filename, 'w', newline='') as file:
        writer = csv.writer(file, delimiter=' ')
        writer.writerow(["N=%d"%N])
        writer.writerow(["err:", err])
        writer.writerow(["finDiff", "analytic"])
        
        for i in range(length):
            writer.writerow([fd_ps[i], al_ps[i]])
        
        print("  saved csv")
        file.close()


# filename = "%d_linear_pressure_%d"%(n, N)

# write_solution(filename, domain.Nx, fd_pressure.ps, al_pressure.ps, N, infNorm_err)

