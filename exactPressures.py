#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 25 18:24:34 2023

@author: sarahdennis
"""
import domain as dfd
import _graphics as graph
import numpy as np

    
class Pressure:
    
    def __init__(self, domain, ps, p0, pN, p_str):
        self.ps = ps
        self.pxs = dfd.center_diff(ps, domain)
        self.pxxs = dfd.center_second_diff(ps, domain)
        self.p_str = p_str
        self.p0 = p0
        self.pf = pN
            
    def plot(self, domain):
        graph.plot_2D(self.ps, domain.xs, "Exact Pressure (%s)"%self.p_str, "p" )
        
    def plot_all(self, domain):
        graph.plot_2D_multi([self.ps, self.pxs, self.pxxs], domain.xs, "Exact Pressure (%s)"%self.p_str, ["p", "px", "pxx"])

class CorrugatedPressure(Pressure):
    def __init__(self, domain, height, p0, pN):
        p_str = "Corrugated"
        ps = np.zeros(domain.Nx)
        for i in range(domain.Nx):
            h = height.hs[i]
            hx = height.hxs[i]
            #TODO how is this derived?
            ps[i] = -6*domain.eta*domain.U * (h + height.h_mid)/((height.k*height.h_mid)**2*(2 + height.r**2)) * hx / h**2
        super().__init__(domain, ps,p0, pN,  p_str)
        
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

        #---------------
        M = np.zeros((2*n + 1, 2*n + 1))
        
        #B:= top left corner of M, 1/dx diagonal
        B = np.zeros((n+1, n))
        B_diag_neg = [-1/height.step_width]*n
        B_diag_pos = [1/height.step_width]*n
        B[0:n , 0:n] += np.diagflat(B_diag_neg)
        B[1:n+1, 0:n] += np.diagflat(B_diag_pos)  
        
        #C:= bottom right corner of M, hj-hi diagonal
        C = np.zeros((n, n+1))
        C_diag_neg = [-h**3 for h in height.h_steps[0:n]] #h_steps = [h1, h2,...hn+1]
        C_diag_pos = [h**3 for h in height.h_steps[1:n+1]]
        C[0:n, 0:n] += np.diagflat(C_diag_neg)
        C[0:n, 1:n+1] += np.diagflat(C_diag_pos)

        #---------------
        M[0:n+1, 0:n+1] = np.identity(n+1)
        M[0:n+1, n+1:2*n+1] = B
        M[n+1:2*n+1, 0:n+1] = C
        #print("M: \n", M)
        
        #---------------

        rhs = np.zeros(2*n + 1)
        
        rhs[0] = -p0/height.step_width
        rhs[n] = pN/height.step_width
        
        for k in range(n):
            rhs[n+1 + k] = (height.h_steps[k+1] - height.h_steps[k]) * 6*domain.eta*domain.U
        
        #print("rhs: ", rhs)
        
        #---------------
        sol = np.linalg.solve(M, rhs)
        
        p_slopes = sol[0:n+1]
        p_extrema =  sol[n+1:2*n+1]
        #print("slopes: ", p_slopes)
        #print("pressures: ", p_extrema)
        #---------------
        
        #ps:= piecewise-linear p evaluated on grid                                                        
        ps = np.zeros(domain.Nx)
        Lx = int(height.step_width/domain.dx) # num xs per step 

        xa = domain.x0 
        xb = domain.x0 + height.step_width
        pa = p0
        for i in range(n+1):
            ps[i*Lx:(i+1)*Lx] = line(xa, xb, Lx, domain.dx, p_slopes[i], pa)
            xa = xb
            xb += height.step_width
            if i != n:
                pa = p_extrema[i]
        super().__init__(domain, ps, p0, pN, p_str)

    
def line(x0, xf, Lx, dx, m, y0): 
    ys = np.zeros(Lx)
    for i in range(Lx): 
        xi = x0 + i * dx
        ys[i] = m*(xi-x0) + y0
    return ys # [y0, ..., yf]
        

        
