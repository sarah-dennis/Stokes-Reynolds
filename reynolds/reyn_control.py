        #!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 22 15:05:40 2023
@author: sarahdennis
"""
import reyn_examples as exs
import numpy as np

# -- solver called on initialization: self.ps = *.solve()

example = exs.Random()
# example = exs.Constant_Ex1()
# example = exs.Linear_Ex1()       
# example = exs.Step_Ex1()
# example = exs.PiecewiseLinear_Ex1()
# example = exs.PiecewiseLinear_Ex2()

fd_solution = exs.Custom(example) # use hs from another example
#------------------------------------------------------------------------------
# plotting 
#------------------------------------------------------------------------------

example.plot_p()
# fd_solution.plot_p()
example.plot_v()

#------------------------------------------------------------------------------
# Error
#------------------------------------------------------------------------------

infNorm_err = np.max(np.abs(example.ps - fd_solution.ps))
print("Analytic to Numerical Error: %.8f"%infNorm_err)

#------------------------------------------------------------------------------
# CSV writing
#------------------------------------------------------------------------------

# def write_solution(filename, length, fd_ps, al_ps, err):

#     with open(filename, 'w', newline='') as file:
#         writer = csv.writer(file, delimiter=' ') 
#         writer.writerow(["err:", err])
#         writer.writerow(["finDiff", "analytic"])
        
#         for i in range(length):
#             writer.writerow([fd_ps[i], al_ps[i]])
        
#         print("  saved csv")
#         file.close()

# dm = example.solver.height
# filename = "%s_pwLinear_pressure_%d"%(example.solver.filename,dm.Nx )

# write_solution(filename, dm.Nx, fd_solution.ps, example.ps, infNorm_err)

