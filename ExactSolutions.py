#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 25 18:24:34 2023

@author: sarahdennis
"""
import domain as dfd
import _graphics as graph
import numpy as np
    
class ExactPressure:
    
    def __init__(self, domain, ps, p_str):
        self.ps = ps
        self.pxs = dfd.center_diff(ps, domain)
        self.pxxs = dfd.center_second_diff(ps, domain)
        self.p_str = p_str
        self.p0 = ps[0]
        self.pf = ps[-1]
            
    def plot(self, domain):
        graph.plot_2D(self.ps, domain.xs, "Exact Pressure (%s)"%self.p_str, "p" )
        
    def plot_all(self, domain):
        graph.plot_2D_multi([self.ps, self.pxs, self.pxxs], domain.xs, "Exact Pressure (%s)"%self.p_str, ["p", "px", "pxx"])


class CorrugatedPressure(ExactPressure):
    def __init__(self, domain, height, p0, pN):
        p_str = "Corrugated"
        ps = np.zeros(domain.Nx)
        for i in range(domain.Nx):
            h = height.hs[i]
            hx = height.hxs[i]
            #TODO how is this derived?
            ps[i] = -6*domain.eta*domain.U * (h + height.h_mid)/((height.k*height.h_mid)**2*(2 + height.r**2)) * hx / h**2
        super().__init__(domain, ps, p_str)
        
class WedgePressure(ExactPressure):
    def __init__(self, domain, height, p0, pN):
        p_str = "Wedge"
        ps = np.zeros(domain.Nx)
        a = height.h_max/height.h_min
        L = domain.xf - domain.x0
        
        for i in range(domain.Nx):
            X = domain.xs[i]/L
            Pi = a/(1-a**2)*(1/self.H(X, a)**2 - 1/a**2) - 1/(1-a)*(1/self.H(X, a)-1/a)
            ps[i] = Pi * 6 * domain.eta * domain.U * L / height.h_min**2
        super().__init__(domain, ps, p_str)
    def H(self, X, a):
        return a + (1-a)*X

# class StepPressure(ExactPressure):

#     def __init__(self, domain, height, p0, pN):
#         p_str = "Step"

#         ps = np.zeros(domain.Nx)
        
#         for i in range(domain.Nx):

#             hl = height.h_left
#             hr = height.h_right
#             ll = height.l_left
#             lr = height.l_right
            
#             m_in = 6*domain.eta*domain.U* (hl - hr) / (hl**3 + hr**3 * ll/lr) 
#             m_out = (-ll / lr) * m_in
            
            
#             if domain.xs[i] <= height.x1:
#                 ps[i] = m_in * (domain.xs[i]-domain.xs[0]) + p0
#             else:
#                 ps[i] = m_out *(domain.xs[i]-domain.xs[-1]) + pN
        
#         super().__init__(domain, ps, p_str)
        

class TwoStepPressure(ExactPressure):

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
        
        super().__init__(domain, ps, p_str)



class SquareWavePressure(ExactPressure):
    
    def __init__(self, domain, height, p0, pN, L):
        #period L
        
        p_str = "%d-Step"%height.n_step

        rhs_ps = np.zeros(height.n_step)
        rhs_ps[0] = -p0/domain.dx
        rhs_ps[height.n_step-1] = pN/domain.dx
        
        rhs_ss = np.zeros(n_step - 1)
        
        for i in range(n_step-1):
 
            
            rhs_ss[i] = 6*domain.eta*domain.U * (h2 - h1)
        
        A = np.identity(n_step)
        D = np.zeros((n_step-1, n_step-1))
        
        B = np.zeros((n_step, n_step-1))
        B_diag_u = [-1/domain.dx]*(n_step-1) 
        B_diag_l = [1/domain.dx]*(n_step-1)
        B[:n_step-1][:]+= np.diagflat(B_diag_u)
        B[1:][:] += np.diagflat(B_diag_l)  
        
    
        C = np.zeros((n_step-1, n_step))
        C_diag = [h**3 for h in height.hs]
        C_diag_u = [-h**3 for h in height.hs[1:]]
        C += np.diagflat(C_diag)
        C[:][1:]+= np.diagflat(C_diag_u)
        C[-1][0] = -height.hs[0]**3

        M = np.zeros((n_step * (n_step-1), n_step * (n_step-1)))
        M[0:n_step][0:n_step] = A
        M[0:n_step][n_step:] = B
        M[n_step:][0:n_step] = C
        M[n_step:][n_step:] = D
        

        
        
        rhs = np.zeros(2*n_step - 1)
        rhs[0:n_step] = rhs_ss
        rhs[n_step:] = rhs_ps
        
        sol = np.linalg.solve(M, rhs)
        
        super().__init__(domain, sol, p_str)

        
        





















        
        
