# -*- coding: utf-8 -*-
"""
Created on Fri May 24 08:36:17 2024

@author: sarah
"""
import numpy as np
import graphics
import stokes_readwrite as rw
from scipy.signal import argrelextrema as relEx
#-------------------------------------------------------------------------------

def compare_Ns(ex, N_min, Ns, N_max):
    ex_min = ex(N_min)

    ex_max = ex(N_max)
    u_max, v_max, psi_max, past_iters = rw.read_solution(ex_max.filestr+".csv", ex_max.Nx * ex_max.Ny)
    psi_max = psi_max.reshape((ex_max.Ny,ex_max.Nx))
    mult_max = int(N_max/N_min)
    
    err = np.zeros((len(Ns)+1,ex_min.Ny,ex_min.Nx))
    err_inf = np.zeros(len(Ns)+1)
    err_l1 = np.zeros(len(Ns)+1)
    err_l2 = np.zeros(len(Ns)+1)
    
    for n in range(len(Ns)+1):
        if n == 0:
            ex_n = ex(N_min)
            u, v, psi_n, past_iters = rw.read_solution(ex_n.filestr+".csv", ex_n.Nx * ex_n.Ny)
            psi_n=psi_n.reshape((ex_min.Ny, ex_min.Nx))
            mult = 1
        else:
            ex_n = ex(Ns[n-1])
            u, v, psi_n, past_iters = rw.read_solution(ex_n.filestr+".csv", ex_n.Nx * ex_n.Ny)
            psi_n=psi_n.reshape((ex_n.Ny, ex_n.Nx))
            mult = int(Ns[n-1]/N_min)
        
        for k_min in range(ex_min.Ny*ex_min.Nx):
        # all indices (i,j) on grid N_min
            i = k_min % ex_min.Nx
            j = (k_min // ex_min.Nx)
            
            # -> indices on grid N || N_min
            i_n = mult * i
            j_n = mult * j

            # -> indices on grid N_max || N_min
            i_max = mult_max * i
            j_max = mult_max * j
            
            
            norm =  abs(psi_max[j_max,i_max])
            err[n,j,i] = abs(psi_max[j_max,i_max] - psi_n[j_n,i_n])
            
            if not np.isclose(norm,0) :
                err[n,j,i] /= norm

            err_l1[n] += err[n,j,i]
            err_l2[n] += err[n,j,i]**2
        # graphics.plot_contour_mesh(err[n], ex_min.xs, ex_min.ys, "err", ['x', 'y', 'err'], True)   
        err_inf[n] = np.max(err[n])
        err_l2[n] = np.sqrt(err_l2[n])
        
    l1_rate = convg_rate(err_l1)
    l2_rate = convg_rate(err_l2)
    inf_rate = convg_rate(err_inf)

    print(f"{np.array2string(l1_rate, precision=2)}")
    print(f"{np.array2string(l2_rate, precision=2)}")
    print(f"{np.array2string(inf_rate, precision=2)}")

    return err_l1, err_l2, err_inf
        
def convg_rate(errs, order=2):
    n = len(errs)
    rates = np.zeros(n-1)
    for k in range(n-1):
        rates[k]=errs[k+1]/(errs[k]**order)
    return rates
#------------------------------------------------------------------------------
# for triangle examples:
#------------------------------------------------------------------------------
def get_criticals(tri, N):
    tri = tri(N)
    u, v, psi, past_iters = rw.read_solution(tri.filestr+".csv", tri.Nx * tri.Ny)
    
    left, right = get_boundary(tri, psi)
    center = get_center(tri, psi)
    
    left_signs = []
    right_signs = []
    
    center_extrs = []

    extr_inds = relEx(abs(center[:,2]), np.greater)[0]
  
    for i in extr_inds:
        x,y,p = center[i]
        # print("(x:%.1f, y:%.6f) p=%.5e"% (x,y,p))
        center_extrs.append([x,y,p])
    
    
    sign_ref = 0
    for (x, y, p) in left:
        sign_new = np.sign(p)
        if sign_new != 0 and sign_new != sign_ref:
            sign_ref = sign_new
            left_signs.append([x,y])

    sign_ref = 0
    for (x, y, p) in right:
        sign_new = np.sign(p)
        if sign_new != 0 and sign_new != sign_ref:
            sign_ref = sign_new
            right_signs.append([x,y])
    
    # order starting at y=yL not y=y0
    left_signs.reverse()
    right_signs.reverse()
    center_extrs.reverse()
    
    return center_extrs, left_signs, right_signs

def get_boundary(tri, psi):
    n = tri.Nx
    m = tri.Ny
    
    i_mid = n//2
    
    left = np.zeros((m,3)) # h, (x, y, psi)
    right = np.zeros((m,3)) 
    
    
    for j in range(m):
        
        y = tri.y0 + j*tri.dx
        
        dj = j % tri.slopes[0]
        di = int(j//tri.slopes[0]) 
        
        if dj == 0: #true boundary points
            k_left = j*n + i_mid - di 
            k_right = j*n + i_mid + di
        else: #interior boundary points
            k_left = j*n + i_mid - di + 1  
            k_right = j*n + i_mid + di - 1
        
        x_left = (i_mid - di + 1)*tri.dx + tri.x0
        x_right = tri.xf - (i_mid - di - 1)*tri.dx
        
        left[j] = [x_left, y, psi[k_left]]
        right[j] = [x_right, y, psi[k_right]]
        
    return left, right

def get_center(tri, psi):
    n = tri.Nx 
    m = tri.Ny 
    
    i_mid = n//2
    x = i_mid * tri.dx+ tri.x0
    
    center = np.zeros((m,3)) # x, y, psi
    
    for j in range(m):
        
        y = tri.y0 + j*tri.dx
        
        k = j*n + i_mid
        
        center[j] = [x, y, psi[k]]
        
    return center


#------------------------------------------------------------------------------

def tri_compare_N(tri, Ns, N_max): #Ns: [44, 120, 240, 512, 1000]
    
    tru_extrs, tru_left, tru_right = get_criticals(tri, N_max)
    
    M = len(Ns)
    a = len(tru_extrs)
    c = len(tru_left)
    d = len(tru_right)

    err_extrs = np.zeros((M, a)) # error in stream extrema along x=0.5
    err_extrs_y = np.zeros((M, a)) # error in stream extrema along x=0.5
    err_left = np.zeros((M, c))  # error in saddle point y along left boundary
    err_right = np.zeros((M, d)) #  ''        ''         ''      right   ''
    
    for i in range(M):
        N = Ns[i]
        
        N_extrs,  N_left, N_right = get_criticals(tri, N)

        for j in range(a):
            if j < len(N_extrs): #[0:x, 1:y, 2:psi]                 
                err_extrs[i, j] = np.abs(tru_extrs[j][2] - N_extrs[j][2])
            else:
                err_extrs[i, j] = None
                
        for j in range(a):
            if j < len(N_extrs): #[0:x, 1:y, 2:psi]                 
                err_extrs_y[i, j] = np.abs(tru_extrs[j][1] - N_extrs[j][1])
            else:
                err_extrs_y[i, j] = None
                
        for j in range(c):
            if j < len(N_left): #[0:x, 1:y]
                err_left[i, j] = np.abs(tru_left[j][1] - N_left[j][1])
            else:
                err_left[i, j] = None
        
        for j in range(d):
            if j < len(N_right):#[0:x, 1:y]
                err_right[i, j] = np.abs(tru_right[j][1] - N_right[j][1])
            else:
                err_right[i, j] = None 
        
    return err_extrs.T, err_extrs_y.T, err_left.T, err_right.T
   
#------------------------------------------------------------------------------ 
# For BFS examples
#------------------------------------------------------------------------------
def get_attatchments(ex, N):

    step = ex(N)
    u, v, psi, past_iters = rw.read_solution(step.filestr+".csv", step.Nx * step.Ny)
    
    psi_2D = psi.reshape((step.Ny,step.Nx))
    psi_xs_yf = psi_2D[step.jf_out] #reattatchment wall 
    psi_ys_xstep = psi_2D[:,step.i_step] # detatchment wall
    
    xs_saddle = []
    sign_ref = 0
    for i in range(1, step.Nx):
        sign_new = np.sign(psi_xs_yf[i])
        if sign_new != sign_ref:
            sign_ref = sign_new
            xs_saddle.append(step.xs[i])

    
    ys_saddle = []
    sign_ref = 0
    for j in range(1, step.Ny):
        sign_new = np.sign(psi_ys_xstep[j])
        if sign_new != sign_ref:
            sign_ref = sign_new
            ys_saddle.append(step.ys[j])

        
    return xs_saddle, ys_saddle
    

#------------------------------------------------------------------------------

def bfs_compare_N(ex, Ns, N_max):
    
    tru_xs, tru_ys = get_attatchments(ex, N_max)
    
    M = len(Ns)
    a = len(tru_xs)
    b = len(tru_ys)

    err_xs = np.zeros((M, a)) # error in stream extrema along y=0
    err_ys = np.zeros((M, b)) # error in stream extrema along x=x_step
    
    for i in range(M):
        N = Ns[i]
        
        N_xs, N_ys = get_attatchments(ex, N)

        for j in range(a):
            if j < len(N_xs): #[0:x, 1:y, 2:psi]                 
                err_xs[i, j] = np.abs(tru_xs[j] - N_xs[j])
            else:
                err_xs[i, j] = None
                
        for j in range(b):
            if j < len(N_ys): #[0:x, 1:y, 2:psi]                 
                err_ys[i, j] = np.abs(tru_ys[j] - N_ys[j])
            else:
                err_ys[i, j] = None
 
        
    return err_xs.T, err_ys.T




