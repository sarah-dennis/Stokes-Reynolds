#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 25 18:24:34 2023

@author: sarahdennis
"""
import Reynolds_1D as ry
import _graphics as graph
import numpy as np
import schur
import time
import squareWave as sw

class Pressure:
    
    def __init__(self, domain, ps, p0, pN, p_str, time=0):
        self.ps = ps
        #self.pxs = dfd.center_diff(ps, domain)
        #self.pxxs = dfd.center_second_diff(ps, domain)
        self.p_str = p_str
        self.p0 = p0
        self.pf = pN
        self.time = time
            
    def plot(self, domain):
        graph.plot_2D(self.ps, domain.xs, "Analytic Pressure (%s)"%self.p_str, "p", "x")
        
        
class ReynoldsPressure(Pressure):
    def __init__(self, domain, height, p0, pN):
        p_str = "Reynolds"
        ps = ry.solve(domain, height, p0, pN)
        super().__init__(domain, ps, p0, pN, p_str)
        

class CorrugatedPressure(Pressure): #sinusoidal
    def __init__(self, domain, height, p0, pN):
        p_str = "Corrugated"
        ps = np.zeros(domain.Nx)
        for i in range(domain.Nx):
            h = height.hs[i]
            hx = height.hxs[i]
            #TODO how is this derived?
            ps[i] = -6*domain.eta*domain.U * (h + height.h_mid)/((height.k*height.h_mid)**2*(2 + height.r**2)) * hx / h**2
        super().__init__(domain, ps, p0, pN,  p_str)
        
class WedgePressure(Pressure):
    def __init__(self, domain, height, p0, pN):
        p_str = "Wedge"
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
        p_str = "Step"

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
        

class TwoStepPressure(Pressure):

    def __init__(self, domain, height, p0, pN):
        p_str = "Two Step"

        ps = np.zeros(domain.Nx)
        
        h1 = height.h1
        h2 = height.h2
        h3 = height.h3
        l1 = height.l1
        l2 = height.l2
        l3 = height.l3
        
        #p2 is pressure at end of second step
        p2_a = (h1**3 * l2/l1 + h2**3) * (h3**3 * pN/l3 + 6*domain.U*domain.eta * (h2-h3))
        p2_b = (h2**3 + h3**3 * l2/l3) * (h1**3 * p0/l1 - 6*domain.U*domain.eta * (h2-h1))
        p2_c = h3**3/l3 * (h1**3 * l2/l1 + h2**3) * (p0 + h3**3/h1**3 * l1/l3 *pN - 6*domain.U*domain.eta * l1/h1**3 * (h3-h1))
        p2_d = h1**3/l1 * (h2**3 + h3**3 * l2/l3) - h3**6/l3**2 * l1/h1**3 * (h1**3 * l2/l1 + h2**3)
        p2 = (p2_a + p2_b - p2_c)/p2_d

        #p1 is pressure at end of first step
        p1 = p0 + h3**3 / h1**3 * l1/l3 * (pN-p2) - 6*domain.U*domain.eta * l1/h1**3 * (h3-h1)

        #pressure slopes
        m1 = (p1-p0)/l1
        m2 = (p2-p1)/l2
        m3 = (pN-p2)/l3 
        
        for i in range(domain.Nx):

            if domain.xs[i] <= height.x1:
                ps[i] = m1 * (domain.xs[i]-domain.x0) + p0
                
            elif domain.xs[i] <= height.x2:
                ps[i] = m2 *(domain.xs[i]-height.x1) + p1
                
            else:
                ps[i] = m3 *(domain.xs[i]-height.x2) + p2
        
        super().__init__(domain, ps, p0, pN, p_str)



class SquareWavePressure(Pressure):

    def __init__(self, domain, height, p0, pN):
        n = height.n_steps
        p_str = "%d-Step Square Wave"%n
        rhs = sw.make_RHS(domain, height, p0, pN)
        
    # 1. Build M, python-solve M @ sol = rhs
        t0 = time.time()
        M = sw.make_M(domain, height, p0, pN)     
        t1 = time.time()
        
        sol_1 = np.linalg.solve(M, rhs)
        t2 = time.time()
        
        print("Time building M: %.5f "%(t1-t0))
        print("Time solving M x = b : %.5f "%(t2-t1))
        print("Total time with M: %.5f  \n"%(t2-t0))

    # 2. Build S, evaluate M_inv(S) @ rhs = sol 
        t4 = time.time()
        K_off_diag, K_center_diag = sw.make_schurCompDiags(height)
        phis = schur.get_phis(n, K_off_diag, K_center_diag)
        thetas = schur.get_thetas(n, K_off_diag, K_center_diag)
        K_off_diag_prod = schur.triDiagProd(K_off_diag)
        t5 = time.time()
        
        sol_2 = np.zeros(2*n+1)
        for i in range(2*n+1):
            sol_2[i] = sw.lhs_i(rhs, height, thetas, phis, K_off_diag_prod, i)
        t9 = time.time()
        
        print("Time building K, phis, thetas: %.5f"%(t5 - t4))
        print("Time evaluating M^-1 @ b = x: %.5f"%(t9 - t5))
        print("Total time with M^-1: %.5f \n"%(t9 - t4))
        
        #---------------
        p_slopes = sol_2[0:n+1]
        p_extrema = sol_2[n+1:2*n+1]

        ps = self.make_ps(domain, height, p0, pN, p_slopes, p_extrema)
        
        eval_time = t9-t4
        super().__init__(domain, ps, p0, pN, p_str, eval_time)
 
    #Construct piecewise linear pressure on Nx grid
    def make_ps(self, domain, height, p0, pN, slopes, extrema):
        ps = np.zeros(domain.Nx)
        L = height.step_width
        
        k = 0
        x_k = domain.x0
        p_k = p0
        slope_k = slopes[k]

        for i in range(domain.Nx-1):
            x = domain.xs[i]
            
            #if x is in a new step
            if x > domain.x0 + (k+1)*L:
                k += 1
                x_k = domain.x0 + k*L
                p_k = extrema[k-1]
                slope_k = slopes[k]

            ps[i] = slope_k*(x-x_k) + p_k
        
        ps[-1] = pN
        return ps

    
        
        
        
        
        
        
        

