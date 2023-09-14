#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 25 18:24:34 2023

@author: sarahdennis
"""

import _graphics as graph
import numpy as np
import schur
import time
import squareWave as sw
import domain as dfd

from scipy.sparse.linalg import gmres

class Pressure:
    
    def __init__(self, domain, ps, p0, pN, p_str, time=0):
        self.ps = ps
        self.pxs = dfd.center_diff(self.ps, domain)
        self.pxxs = dfd.center_second_diff(self.ps, domain)
        self.p_str = p_str
        self.p0 = p0
        self.pf = pN
        self.time = time
            
    def plot(self, domain):
        graph.plot_2D(self.ps, domain.xs, "Pressure (%s)"%self.p_str, "Pressure $p(x)$", "$x$")
        
class ConstantPressure(Pressure):
    def __init__(self, domain, height, p0, pN):
        p_str = "Analytic"
        ps = [p0 for i in range(domain.Nx)]
        super().__init__(domain, ps, p0, pN,  p_str)
        

class CorrugatedPressure(Pressure): #sinusoidal
    def __init__(self, domain, height, p0, pN):
        p_str = "Analytic"
        ps = np.zeros(domain.Nx)
        for i in range(domain.Nx):
            h = height.hs[i]
            hx = height.hxs[i]
            ps[i] = -6*domain.eta*domain.U * (h + height.h_mid)/((height.k*height.h_mid)**2*(2 + height.r**2)) * hx / h**2
        super().__init__(domain, ps, p0, pN,  p_str)
        
class WedgePressure(Pressure):
    def __init__(self, domain, height, p0, pN):
        p_str = "Analytic"
        ps = np.zeros(domain.Nx)
        a = height.h_max/height.h_min
        L = domain.xf - domain.x0
        for i in range(domain.Nx):
            X = domain.xs[i]/L
            Pi = a/(1-a**2)*(1/self.H(X, a)**2 - 1/a**2) - 1/(1-a)*(1/self.H(X, a)-1/a)
            ps[i] = Pi * 6 * domain.eta * domain.U * L / height.h_min**2
        super().__init__(domain, ps, p0, pN, p_str)
        
    def H(self, X, a):
        return a + (1-a)*X

class StepPressure(Pressure):

    def __init__(self, domain, height, p0, pN):
        p_str = "Analytic"

        ps = np.zeros(domain.Nx)
        
        for i in range(domain.Nx):

            hl = height.h_left
            hr = height.h_right
            ll = height.l_left
            lr = height.l_right
            
            m_in = 6*domain.eta*domain.U* (hl - hr) / (hl**3 + hr**3 * ll/lr) 
            m_out = (-ll / lr) * m_in
            
            
            if domain.xs[i] <= height.x1:
                ps[i] = m_in * (domain.xs[i]-domain.xs[0]) + p0
            else:
                ps[i] = m_out *(domain.xs[i]-domain.xs[-1]) + pN
        
        super().__init__(domain, ps, p0, pN, p_str)
        

class SquareWavePressure_pySolve(Pressure):

    def __init__(self, domain, height, p0, pN):
        n = height.n_steps
        p_str = "numpy-linalg"
        rhs = sw.make_RHS(domain, height, p0, pN)
        
        t0 = time.time()
        M = sw.make_M(domain, height, p0, pN)  

        t1 = time.time()
        
        sol = np.linalg.solve(M, rhs)
        t2 = time.time()
        
        print("\n Python Inverse Solve")
        print("Prep time: %.5f "%(t1-t0))
        print("Solve time: %.5f "%(t2-t1))
        print("Total time: %.5f "%(t2-t0))
        
        condNum = np.linalg.cond(M)
        print("condNum: %.5f"%condNum)
        
        p_slopes = sol[0:n+1]
        p_extrema = sol[n+1:2*n+1]

        ps = sw.make_ps(domain, height, p0, pN, p_slopes, p_extrema)
        
        super().__init__(domain, ps, p0, pN, p_str, t2-t1)
        
class SquareWavePressure_gmresSolve(Pressure):

# TODO: look in to preconditioners
    def __init__(self, domain, height, p0, pN):
        n = height.n_steps
        p_str = "gmres"
        rhs = sw.make_RHS(domain, height, p0, pN)
        
        t0 = time.time()
        
        # M_linOp = sw.swLinOp(n, height.step_width, height.h_steps)
        
        M_linOp = sw.swLinOp_deClass(n, height.step_width, height.h_steps)
        
        t1 = time.time()
        
        tol = 1e-12
        # max_iter = 200
        
        sol, exit_code = gmres(M_linOp, rhs, tol=tol)
    
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
        p_str = "schur-inv"
        rhs = sw.make_RHS(domain, height, p0, pN)

       # Build S, evaluate M_inv(S) @ rhs = sol 
        t0 = time.time()
        center_diag, off_diag = sw.make_schurCompDiags(height)
        C = schur.get_Cs(n, center_diag, off_diag)
        C_prod = schur.triDiagProd(C)
        D = schur.get_Ds(n, center_diag, off_diag, C)
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
 
        
 
class SquareWavePressure_schurLUSolve(Pressure):
    def __init__(self, domain, height, p0, pN):
        n = height.n_steps
        L = height.step_width
        p_str = "LU"
        
        t0 = time.time()
        rhs = sw.make_RHS(domain, height, p0, pN)
        center_diag, off_diag = sw.make_schurCompDiags(height)
        C = schur.get_Cs(n, center_diag, off_diag)
        C_prod = schur.triDiagProd(C)
        D = schur.get_Ds(n, center_diag, off_diag, C)

        t1 = time.time()
        
        # L block -  fwd sub
        # | I  0 | |v| = |f|
        # | B2 K | |w|   |g|
        # {v = f, w = S ( g - B2 f)}
        
        w = np.zeros(n)
        for i in range(n):
            w_i = 0
            for j in range(n):
                s_ij = schur.S_ij(n, C_prod, D, i, j)
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
        center_diag, off_diag = sw.make_schurCompDiags(height)
        C = schur.get_Cs(n, center_diag, off_diag)
        D = schur.get_Ds(n, center_diag, off_diag, C)
        t1 = time.time()
        
        # L block -  fwd sub
        # | I  0 | |v| = |f|
        # | B2 K | |w|   |g|
        # {v = f, w = S ( g - B2 f)}
        
        w = np.zeros(n)
        for i in range(n):
            w_i = 0
            for j in range(n):
                s_ij = schur.S_ij_flops(n, C, D, i, j)
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

