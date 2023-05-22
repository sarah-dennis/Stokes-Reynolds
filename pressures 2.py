#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 25 18:24:34 2023

@author: sarahdennis
"""
import _graphics as graph
import numpy as np
import schur

class Pressure:
    
    def __init__(self, domain, ps, p0, pN, p_str):
        self.ps = ps
        #self.pxs = dfd.center_diff(ps, domain)
        #self.pxxs = dfd.center_second_diff(ps, domain)
        self.p_str = p_str
        self.p0 = p0
        self.pf = pN
            
    def plot(self, domain):
        graph.plot_2D(self.ps, domain.xs, "Analytic Pressure (%s)"%self.p_str, "p", "x")
        
    def plot_diffs(self, domain):
        graph.plot_2D_multi([self.ps, self.pxs, self.pxxs], domain.xs, "Analytic Pressure Derivatives (%s)"%self.p_str, ["p", "px", "pxx"])

class UnknownPressure(Pressure):
    def __init__(self, domain, height, p0, pN):
        p_str = "Exact Pressure Unknown"
        ps = np.zeros(domain.Nx)
        super().__init__(domain, ps, p0, pN,  p_str)

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
    
    # def build_sqrWave_Matrix(self, domain, height, p0, pN):

    #     n = height.n_steps
    #     hs = height.h_steps 
    #     #---------------
    #     M = np.zeros((2*n + 1, 2*n + 1))

    #     #B:= top left corner of M, 1/L diagonal
    #     B = np.zeros((n+1, n))
    #     B_diag_neg = [-1/height.step_width]*n
    #     B_diag_pos = [1/height.step_width]*n
    #     B[0:n , 0:n] += np.diagflat(B_diag_neg)
    #     B[1:n+1, 0:n] += np.diagflat(B_diag_pos)  
        
    #     #C:= bottom right corner of M, hj-hi diagonal
    #     C = np.zeros((n, n+1))
    #     C_diag_neg = [-h**3 for h in hs[0:n]] 
    #     C_diag_pos = [h**3 for h in hs[1:n+1]]
    #     C[0:n, 0:n] += np.diagflat(C_diag_neg)
    #     C[0:n, 1:n+1] += np.diagflat(C_diag_pos)

    #     #...
    #     M[0:n+1, 0:n+1] = np.identity(n+1)
    #     M[0:n+1, n+1:2*n+1] = B
    #     M[n+1:2*n+1, 0:n+1] = C
        
    #     #---------------
    #     rhs = np.zeros(2*n + 1)
        
    #     rhs[0] = -p0/height.step_width
    #     rhs[n] = pN/height.step_width
        
    #     for k in range(n):
    #         rhs[n+1 + k] = (hs[k+1] - hs[k]) * 6*domain.eta*domain.U
        
    #     return M, rhs
    
    def make_Minv_schurComp(self, height, S):
        
        n = height.n_steps
        
        block_11 = np.zeros((n+1, n+1))
        block_12 = np.zeros((n+1, n))
        block_21 = np.zeros((n, n+1))

        for i in range(0, n+1):
            for j in range(0, n+1):
                
                block_11[i,j] = self.Id_B1schurCompInvB2_ij(height, S, i, j)
                
                if i < n:
                    block_21[i, j] = self.schurCompInvB2_ij(height, S, i, j)
                
                if j < n:
                    block_12[i,j] = self.B1schurCompInv_ij(height, S, i, j)
                
        return [block_11, block_21], [block_21, S]
    
    
    
    # get diags of schur complement K (symmetric tri-diagonal)
    def make_schurCompDiags(self, height):
        n = height.n_steps
        hs = height.h_steps

        center_diag = np.zeros(n)
        off_diag = np.zeros(n-1)
        
        for i in range (n):
            center_diag[i] = hs[i]**3 + hs[i+1]**3
            if i < n-1:
                off_diag[i] = -hs[i+1]**3
                
        #add a factor of -1/L when done
        return off_diag, center_diag
    
    
    def schurCompInvB2_ij(self, height, S, i, j):
        n = height.n_steps
        if j == 0:
            return -height.h_steps[j]**3 * S[i,j]
        elif j < n:
            return height.h_steps[j]**3 * (S[i, j-1] - S[i,j])
        else:
            return height.h_steps[j]**3 * S[i, j-1]
        
    def B1schurCompInv_ij(self, height, S, i, j):
        L = height.step_width
        n = height.n_steps
        if i == 0:
            return (-1/L) * (-S[i,j])
        elif i < n:
            return (-1/L) * (S[i-1,j] - S[i,j])
        else:
            return (-1/L) * S[i-1,j]
            
    def Id_B1schurCompInvB2_ij(self, height, S, i, j):
        L = height.step_width
        n = height.n_steps
        
        if i == 0 and j == 0:
            return 1 + (1/L) * height.h_steps[j]**3 * S[i,j]
        
        elif i == 0 and j < n :
            return (-1/L) * height.h_steps[j]**3 * (S[i,j-1] - S[i,j])
        
        elif i == 0 and j == n: 
            return (-1/L) * height.h_steps[j]**3 * (S[i,j-1])
        
        elif i == n and j == 0:
            return (-1/L) * height.h_steps[j]**3 * (S[i-1,j])
        
        elif i == n and j < n:
             return (1/L) * height.h_steps[j]**3 * (S[i-1,j-1] - S[i-1,j])
        
        elif i == n and j == n:
            return 1 + (1/L) * height.h_steps[j]**3 * S[i-1,j-1]
        
        elif i < n and j == 0:
            return (-1/L) * height.h_steps[j]**3 * (S[i-1,j] - S[i,j])

        elif i < n and j == n:
            return (1/L) * height.h_steps[j]**3 * (S[i-1,j-1] - S[i,j-1])
        
        else:
            return (i==j) + (1/L) * height.h_steps[j]**3 * (S[i-1,j-1] - S[i-1,j] - S[i, j-1] + S[i,j])    
    

    def __init__(self, domain, height, p0, pN):

        n = height.n_steps
        hs = height.h_steps
        L = height.step_width
        p_str = "%d-Step Square Wave"%n
        
        ##1. Build M and solve
        #M, rhs = self.build_sqrWave_Matrix(domain, height, p0, pN)
        #sol = np.linalg.solve(M, rhs)
        
        #p_slopes = sol[0:height.n_steps+1]
        #p_extrema = sol[height.n_steps+1:2*height.n_steps+1]


        #2. Build M^-1 using schur comp
        K_off_diag, K_center_diag = self.make_schurCompDiags(height)
        S = (-1/L) * schur.make_symTriInv(n, K_off_diag, K_center_diag)
        M_inv = self.make_Minv_schurComp(height, S)
        
        
        print(M_inv[0][0][0][0])
        #fake load with ps=hs
        super().__init__(domain, hs, p0, pN, p_str)
        
        #print("slopes: ", p_slopes)
        #print("pressures: ", p_extrema)
        #---------------
                                                   
        # ps = np.zeros(domain.Nx)
        
        # k = 0
        # x_k = domain.x0
        # p_k = p0
        # slope_k = p_slopes[k]

        # for i in range(domain.Nx-1):
        #     x = domain.xs[i]

        #     if x > domain.x0 + (k+1)*height.step_width:

        #         k += 1
                
        #         #x_k += height.step_width
        #         x_k = domain.x0 + k*height.step_width
                
        #         p_k = p_extrema[k-1] #p_extrema = [p1, ..., pN-1]
         
        #         slope_k = p_slopes[k] #slopes = [s1, ..., sN-1, sN]

        #     ps[i] = slope_k*(x-x_k) + p_k
        
        # ps[-1] = pN
  
        # super().__init__(domain, ps, p0, pN, p_str)
        
        
 #TODO: Dimple pressure
        
        
        
        
        

