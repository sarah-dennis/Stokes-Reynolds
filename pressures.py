#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 25 18:24:34 2023

@author: sarahdennis
"""

import numpy as np
import time

import pressure_sqrWave as sw
import pressure_sawtooth as st
import pressure_finDiff as fd

from scipy.sparse.linalg import gmres


class Pressure:
    def __init__(self, height, p0, pN, p_solver, p_str):
        self.height = height
        self.p0 = p0
        self.pN = pN
        self.solve = p_solver #pressure.solve(height, p0, pN)
        self.p_str = p_str
    

class FinDiffPressure(Pressure):
    def __init__(self, height, p0, pN):
        p_str = "Finite Difference Reynolds ($N_x = %d$)"%height.Nx
        solver = fd.solve
        super().__init__(height, p0, pN, solver, p_str)
        
class LinearPressure(Pressure):
    def __init__(self, height, p0, pN):
        p_str = "Analytic Reynolds"
        p_solver = LinearPressure.solve
        
        super().__init__(height, p0, pN, p_solver, p_str)
    
    def solve(height, p0, pN):
        
        h = height.h0
        
        etaU = 6*height.U*height.visc
        
        cq = h**3 * (pN - p0)/(height.xf - height.x0) - etaU*h
        
        ps = np.zeros(height.Nx)
        
        for i in range(height.Nx):
            dx = height.xs[i] - height.x0
            
            ps[i] = (cq * dx / h**3) + (etaU * dx / h**2) + p0
        return ps
        

class CorrugatedPressure(Pressure): #sinusoidal
    def __init__(self, domain, height, p0, pN):
        p_str = "Reynolds Analytic"
        #TODO: ignores boundary pressures
        ps = np.zeros(domain.Nx)
        for i in range(domain.Nx):
            h = height.hs[i]
            hx = height.hxs[i]
            ps[i] = -6*domain.eta*domain.U * (h + height.h_mid)/((height.k*height.h_mid)**2*(2 + height.r**2)) * hx / h**2
            
        super().__init__(domain, ps, p0, pN,  p_str)
        

class SawtoothPressure(Pressure):
    def __init__(self, domain, height, p0, pN):
        p_str = "Reynolds Piecewise Analytic"
        
        Xs = height.x_peaks #[x0, x1, ..., xN]
        Hs = height.h_peaks #[h0, h1, ..., hN]
        
        slopes = height.slopes #[(i, i+1): i = 0, ... N-1]
        dxs = height.dxs
        N = len(Xs)-1
       

        rhs = st.make_rhs(N, Hs, slopes, dxs, p0, pN, domain.eta*domain.U)
 

        st_linOp = st.stLinOp(N, Hs, slopes, dxs)
     
        tol = 1e-12
        
        cs, exit_code = gmres(st_linOp, rhs, tol=tol)
        ps = st.make_ps(domain, height, cs, N)
        
      
        super().__init__(domain, ps, p0, pN, p_str)
        
        
class RayleighStepPressure(Pressure):

    def __init__(self, domain, height, p0, pN):
        p_str = "Reynolds Analytic"

        ps = np.zeros(domain.Nx)

        h_in = height.hs[0]
        h_out = height.hs[-1]

        x0 = domain.x0
        xm = height.x_step
        xf = domain.xf
        c = 6*domain.eta*domain.U


        m_in_numer = (h_out/h_in)**3 * (p0 - pN)/(xm - xf) - c * (h_out - h_in)/h_in**3
        m_in_denom = 1 - (h_out/h_in)**3 * (xm - x0)/(xm - xf)
        
        m_in = m_in_numer/m_in_denom

        m_out = ((xm - x0)*m_in + (p0 - pN))/(xm - xf)

        for i in range(domain.Nx):
            if domain.xs[i] <= height.x_step:
                ps[i] = m_in * (domain.xs[i] - x0) + p0
            else:
                ps[i] = m_out * (domain.xs[i] - xf) + pN

        self.m_in = m_in
        self.m_out = m_out
        super().__init__(domain, ps, p0, pN, p_str)
        

class SquareWavePressure_pySolve(Pressure):

    def __init__(self, domain, height, p0, pN):
        n = height.n_steps
        p_str = "Reynolds Piecewise Analytic (np.linalg)"
        rhs = sw.make_RHS(domain, height, p0, pN)
        
        t0 = time.time()
        sw_M = sw.make_M(domain, height, p0, pN)  

        t1 = time.time()
        
        sol = np.linalg.solve(sw_M, rhs)
        t2 = time.time()
        
        print("\n Python Inverse Solve")
        print("Prep time: %.5f "%(t1-t0))
        print("Solve time: %.5f "%(t2-t1))
        print("Total time: %.5f "%(t2-t0))
        
        condNum = np.linalg.cond(sw_M)
        print("condNum: %.5f"%condNum)
        
        p_slopes = sol[0:n+1]
        p_extrema = sol[n+1:2*n+1]

        ps = sw.make_ps(domain, height, p0, pN, p_slopes, p_extrema)
        
        super().__init__(domain, ps, p0, pN, p_str, t2-t1)
        
class SquareWavePressure_schurGmresSolve(Pressure):

# TODO: matvec overcalled, look in to preconditioners
    def __init__(self, domain, height, p0, pN):
        n = height.n_steps
        p_str = "Reynolds Piecewise Analytic (schur & gmres)"
        rhs = sw.make_RHS(domain, height, p0, pN)
        
        t0 = time.time()
        
        sw_linOp = sw.swLinOp(n, height.step_width, height.h_steps)
        
        t1 = time.time()
        
        tol = 1e-12
        # max_iter = 200
        
        sol, exit_code = gmres(sw_linOp, rhs, tol=tol)
    
        t2 = time.time()
        
        print("\n Scipy GMRes Solve")
        print("Prep time: %.5f "%(t1-t0))
        print("Solve time: %.5f "%(t2-t1))
        print("Total time: %.5f "%(t2-t0))
        print(exit_code)
        p_slopes = sol[0:n+1]
        p_extrema = sol[n+1:2*n+1]

        ps = sw.make_ps(domain, height, p0, pN, p_slopes, p_extrema)
        
        super().__init__(domain, ps, p0, pN, p_str, t2-t1)

class SquareWavePressure_schurInvSolve(Pressure):
    def __init__(self, domain, height, p0, pN):
        n = height.n_steps
        p_str = "Reynolds Piecewise Analytic (schur inv)"

       # Build S, evaluate M_inv(S) @ rhs = sol 
        t0 = time.time()
        rhs = sw.make_RHS(domain, height, p0, pN)
        center_diag, off_diag = sw.make_schurCompDiags(height)
        C = sw.get_Cs(n, center_diag, off_diag)
        C_prod = sw.triDiagProd(C)
        D = sw.get_Ds(n, center_diag, off_diag, C)
        t1 = time.time()
        
        sol = np.zeros(2*n+1)
        for i in range(2*n+1):
            sol[i] = sw.schurInvSol_i(rhs, height, C_prod, D,  i)
        t2 = time.time()
        
        print("\n Schur Inv. Solve")
        print("Prep time: %.5f"%(t1 - t0))
        print("Solve time: %.5f"%(t2 - t1))
        print("Total time: %.5f "%(t2 - t0))
        
        p_slopes = sol[0:n+1]
        p_extrema = sol[n+1:2*n+1]

        ps = sw.make_ps(domain, height, p0, pN, p_slopes, p_extrema)

        super().__init__(domain, ps, p0, pN, p_str, t2-t1)
 
        
# Current best method ... 
class SquareWavePressure_schurLUSolve(Pressure):
    def __init__(self, domain, height, p0, pN):
        n = height.n_steps
        
        L = height.step_width
        hs = height.h_steps
        p_str = "Reynolds Piecewise Analytic (schur & LU)"
        
        t0 = time.time()
        
        rhs = sw.make_RHS(domain, height, p0, pN)

        center_diag, off_diag = sw.make_schurCompDiags(height)
        
        C = sw.get_Cs(n, center_diag, off_diag)
        C_prod = sw.triDiagProd(C) # used in schur.S_ij
        D = sw.get_Ds(n, center_diag, off_diag, C)

        t1 = time.time()
        
        # L block -  fwd sub
        # | I  0 | |v| = |f|
        # | B2 K | |w|   |g|
        # {v = f, w = S( g - B2 f)}
        
        w = np.zeros(n)
        
        for i in range(n):

            w_ij = 0
            
            for j in range(n):
                
                s_ij = sw.S_ij(n, C_prod, D, i, j)

                if j == 0:
                    w_ij += s_ij * (rhs[n+1+j] - (1/L) * p0 * hs[j]**3)
                    
                elif j == n-1:
                    w_ij += s_ij * (rhs[n+1+j] -  (1/L) * pN * hs[j+1]**3)

                else: 
                    w_ij += s_ij * rhs[n+1+j]
                    
            w[i] = w_ij
            
        # U block  - back sub
        # | I  B1 | |x| = |v|
        # | 0  I  | |y|   |w|
        # {x = v - B1 w, y = w}
        
        x = np.zeros(n+1)
        
        x[0] = 1/L * (w[0] - p0)
        for i in range(1, n):
            x[i] = 1/L * (w[i] - w[i-1])
            
        x[n] = 1/L * (pN - w[n-1])

        t2 = time.time()

        print("\n Schur LU Solve")
        print("Prep time: %.5f"%(t1 - t0))
        print("Solve time: %.5f"%(t2-t1))
        print("Total time: %.5f"%(t2-t0))
        
        ps = sw.make_ps(domain, height, p0, pN, x, w)
        
        super().__init__(domain, ps, p0, pN, p_str, t2-t1)
 


class SquareWavePressure_schurLUSolve_flops(Pressure):
    def __init__(self, domain, height, p0, pN):
        n = height.n_steps
        L = height.step_width
        p_str = "LU-flops"
        
        t0 = time.time()
        rhs = sw.make_RHS(domain, height, p0, pN)
        center_diag, off_diag = sw.make_schurCompDiags(height)\
        # for flops... Store C & D (N+N) instead of C_prod & D (N^2 + N)
        C = sw.get_Cs(n, center_diag, off_diag)
        D = sw.get_Ds(n, center_diag, off_diag, C)
        t1 = time.time()
        
        # L block -  fwd sub
        # | I  0 | |v| = |f|
        # | B2 K | |w|   |g|
        # {v = f, w = S ( g - B2 f)}
        
        w = np.zeros(n)
        for i in range(n):
            w_i = 0
            for j in range(n):
                s_ij = sw.S_ij_flops(n, C, D, i, j)
                w_i += s_ij * (rhs[n+1+j] - 1/L * rhs[j] * height.h_steps[j]**3)
            w[i] = w_i

        # U block  - back sub
        # | I  B1 | |x| = |v|
        # | 0  I  | |y|   |w|
        # {x = v - B1 w, y = w}
        
        x = np.zeros(n+1)
        x[0] = 1/L * (w[0] - p0)
        for i in range(1, n):
            x[i] = 1/L * (w[i] - w[i-1])
        x[n] = 1/L * (pN - w[n-1])
        
        t2 = time.time()

        print("\n Schur LU-flops Solve")
        print("Prep time: %.5f"%(t1 - t0))
        print("Solve time: %.5f"%(t2-t1))
        print("Total time: %.5f"%(t2-t0))
        
        ps = sw.make_ps(domain, height, p0, pN, x, w)
        super().__init__(domain, ps, p0, pN, p_str, t2-t1)

