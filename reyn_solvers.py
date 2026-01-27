# -*- coding: utf-8 -*-
"""
Created on Tue Jan 27 10:10:36 2026

@author: sarah
"""
import numpy as np
from time import time

import domain as dm
import reyn_boundary as bc
import reyn_velocity as rv
import reyn_velocity_ELT as eltv
import reyn_pressure as rp
import reyn_pressure_ELT as eltp
import reyn_perturbed as rpert
from reyn_solution import Reyn_Solution
from reyn_heights import PWC_Height, PWL_Height, make_PWC, make_PWL


class Reynolds_Solver: 
    def __init__(self, Example, BC, args=None):
        self.Example = Example #initialize height = Example(args) in the solver
        self.args = args
        
        self.BC = BC        
        
#----------------------------------------------------------------------------------
    def fd_solve(self, N):
        solver_title = "Reynolds"
        
        height = self.Example(self.args, N)

        t0 = time()
        height.hxs = dm.center_diff(height.hs, height.Nx, height.dx)
        pressure = rp.FinDiff_ReynPressure(height, self.BC)
        tf = time()
            
        velocity = rv.Reyn_Velocity(height, self.BC, pressure)
        t = tf-t0
        
        solution = Reyn_Solution(height, self.BC, pressure, velocity, solver_title, t)

        return solution
#----------------------------------------------------------------------------------
#---------------------------------Piecewise----------------------------------------
#----------------------------------------------------------------------------------
    def pwc_schur_solve(self, N):
        solver_title = "Reynolds"
        
        height = self.Example(self.args, N)

        if not isinstance(height, PWC_Height):
            height = make_PWC(height)
    
        
        t0 = time()
        pressure = rp.PwcSchur_ReynPressure(height, self.BC)
        tf = time()
            
        height.hxs = np.zeros(height.Nx)
        velocity = rv.Reyn_Velocity(height, self.BC,pressure)
        
        t=tf-t0
        solution = Reyn_Solution(height, self.BC, pressure, velocity, solver_title, t)

        return solution
#----------------------------------------------------------------------------------   
   
    def pwc_schur_parallel_solve(self, N):
        solver_title = "Reynolds"
        
        height = self.Example(self.args, N)

        if not isinstance(height, PWC_Height):
            height = make_PWC(height)
        
        t0 = time()
        pressure = rp.PwcSchur_parallel_ReynPressure(height, self.BC)
        tf = time()
        
        velocity = rv.Reyn_Velocity(height, self.BC,pressure)
        

        height.hxs = dm.center_diff(height.hs, height.Nx, height.dx)
        velocity = rv.Reyn_Velocity(height, self.BC, pressure)
        
        t = tf-t0
        solution = Reyn_Solution(height, self.BC, pressure, velocity, solver_title, t)

        return solution

#----------------------------------------------------------------------------------   
    def pwl_schur_solve(self, N):
        solver_title = "Reynolds" #" Piecewise Linear"
        
        height = self.Example(self.args,N)
        
        if not (isinstance(self.BC, bc.Mixed)): #TODO
            raise TypeError('Only prescribed flux for pwl schur solver')
        
        if not (isinstance(height, PWL_Height) or isinstance(height, PWC_Height)): # PWC is a PWL
            height = make_PWL(height)

        t0 = time()
        pressure = rp.PwlSchur_ReynPressure(height, self.BC)
        tf = time()
        t = tf-t0
        height.hxs = dm.center_diff(height.hs, height.Nx, height.dx)
        velocity = rv.Reyn_Velocity(height, self.BC, pressure)
        
        
        solution = Reyn_Solution(height, self.BC, pressure, velocity, solver_title, t)

        return solution
        return pressure, velocity, tf-t0
    
#----------------------------------------------------------------------------------       
    def pwl_gmres_solve(self, N):
        solver_title = "Reynolds" #" Piecewise Linear"
        
        height = self.Example(self.args,N)
        
        t0 = time()
       
        if not (isinstance(self.BC, bc.Mixed)): #TODO
            raise TypeError('Only prescribed flux for pwl gmres solver')
 
       
        if not (isinstance(height, PWL_Height) or isinstance(height, PWC_Height)) :
            height = make_PWL(height)
       
        pressure = rp.PwlGMRes_ReynPressure(height, self.BC)
        tf = time()
        
        # print('pwl gmres time: ', tf-t0)
        height.hxs = dm.center_diff(height.hs, height.Nx, height.dx)
        velocity = rv.Reyn_Velocity(height, self.BC, pressure)
            
        t = tf-t0
        solution = Reyn_Solution(height, self.BC, pressure, velocity, solver_title, t)

        return solution
    
    
#----------------------------------------------------------------------------------
#--------------------ELT--------------------------------------------------------------
#----------------------------------------------------------------------------------

    def fd_TG_ELT_solve(self, N):
        solver_title = "T.G.-ELT"
        t0 = time()
        height = self.Example(self.args, N)
        height.hxs = dm.center_diff(height.hs, height.Nx, height.dx)
        height.h2xs = dm.center_second_diff(height.hs, height.Nx, height.dx)
        height.h3xs = dm.center_third_diff(height.hs, height.Nx, height.dx)
        
        reyn_pressure = rp.FinDiff_ReynPressure(height, self.BC)
        pressure = eltp.TG_ELT_Pressure(height, self.BC, reyn_pressure)

        tf = time()        
        velocity = rv.Reyn_Velocity(height, self.BC, pressure)
        
        t = tf-t0
        solution = Reyn_Solution(height, self.BC, pressure, velocity, solver_title, t)

        return solution
    
    def fd_VA_ELT_solve(self, N):
        solver_title = "VA-ELT"
        t0 = time()
        height = self.Example(self.args, N)
        height.hxs = dm.center_diff(height.hs, height.Nx, height.dx)
        height.h2xs = dm.center_second_diff(height.hs, height.Nx, height.dx)
        height.h3xs = dm.center_third_diff(height.hs, height.Nx, height.dx)
        
        reyn_pressure = rp.FinDiff_ReynPressure(height, self.BC)
        pressure = eltp.VA_ELT_Pressure(height, self.BC, reyn_pressure)
        tf = time()
        velocity = eltv.ELT_Velocity(height,self.BC, pressure)
        
        t = tf-t0
        solution = Reyn_Solution(height, self.BC, pressure, velocity, solver_title, t)

        return solution
    

    def fd_pert_solve(self, N, order, get_both=True):
        height = self.Example(self.args, N)
        
        t0 = time()
        
        height.hxs = dm.center_diff(height.hs, height.Nx, height.dx)
        height.h2xs = dm.center_second_diff(height.hs, height.Nx, height.dx)
        height.h3xs = dm.center_third_diff(height.hs, height.Nx, height.dx)
        
        reyn_pressure = rp.FinDiff_ReynPressure(height, self.BC)
        reyn_velocity = rv.Reyn_Velocity(height, self.BC, reyn_pressure)                   
        t_r = time() -t0
        reyn_sol = Reyn_Solution(height, self.BC, reyn_pressure, reyn_velocity, 'Reynolds', t_r)

        pert = rpert.Perturbed_Solution(height, self.BC, order, reyn_sol)
        
        tf = time()
      
        t = tf-t0
        
        if order >= 2:
            solver_title = "$\epsilon^2$ PLT"
            solution_2 = Reyn_Solution(height, self.BC, pert.pert2_pressure, pert.pert2_velocity, solver_title, t)
            
        if order >= 4:
            solver_title = "$\epsilon^4$ PLT"
            solution_4 = Reyn_Solution(height, self.BC, pert.pert4_pressure, pert.pert4_velocity, solver_title, t)
        
        if order == 2:
            return solution_2
        
        elif order == 4:
            if get_both:
                return solution_2, solution_4
            else:
                return solution_4

        
#----------------------------------------------------------------------------------
    def exact_sinusoid_sol(self, N):
        height = self.Example(self.args, N)
       
        ps = np.zeros(height.Nx)
        t0 = time()
        height.hxs = dm.center_diff(height.hs, height.Nx, height.dx)
        h_recip_dx = [height.h_recip_deriv_fun(x) for x in height.xs]
        for i in range(height.Nx):
            h = height.hs[i] 
            ps[i] = -6*self.BC.U * (h + height.H)/((height.k * height.H)**2 * (2+height.delta**2)) * h_recip_dx[i]
        
        pressure = rp.Reyn_Pressure(height, ps)
        velocity = rv.Reyn_Velocity(height, self.BC, pressure)
        tf = time()
        t = tf-t0
        solver_title = "exact solution"
        solution = Reyn_Solution(height, self.BC, pressure, velocity, solver_title, t)
        return solution