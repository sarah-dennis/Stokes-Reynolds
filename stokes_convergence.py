# -*- coding: utf-8 -*-
"""
Created on Fri May 24 08:36:17 2024

@author: sarah
"""
import numpy as np
import graphics
import stokes_examples as examples
import stokes_readwrite as rw
from scipy.signal import argrelextrema as relEx

#-------------------------------------------------------------------------------
# for triangle examples:
def get_criticals(N):
    tri = examples.biswasEx(N)
    u, v, psi, past_iters = rw.read_solution(tri.filename+".csv", tri.Nx * tri.Ny)
    
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
    
    left = np.zeros((m,3)) # x, y, psi
    right = np.zeros((m,3)) 
    
    
    for j in range(m):
        
        y = tri.y0 + j*tri.dx
        
        dj = j % tri.slope
        di = int(j//tri.slope) 
        
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

# For BFS examples

def get_attatchments(N):
    example = examples.bfs_biswasLowerRe

    step = example(N)
    u, v, psi, past_iters = rw.read_solution(step.filename+".csv", step.Nx * step.Ny)
    
    psi_2D = psi.reshape((step.Ny,step.Nx))
    psi_xs_yf = psi_2D[step.jf_out+1] #reattatchment wall 
    psi_ys_xstep = psi_2D[:,step.i_step+1] # detatchment wall
    
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

def tri_compare_N(Ns, N_max): #Ns: [44, 120, 240, 512, 1000]
    tru_extrs, tru_left, tru_right = get_criticals(N_max)
    
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
        
        N_extrs,  N_left, N_right = get_criticals(N)

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
   
# plot_compare_N([120,240,512,1000],2000)     
def plot_tri_compare_N(Ns, N_max):
    err_extrs, err_extrs_y, err_left, err_right = tri_compare_N(Ns, N_max)
    
    title_extrs_y = "Error to $N^{*}=$%d in $y$ of vortex center"%N_max
    title_extrs = "Error to $N^{*}=$%d in $\psi$ of vortex center"%N_max
    ax_labels_stream = ["N", "$|\psi_{N^{*}} - \psi_{N}|$"]
    ax_labels_y = ["N", "$|y_{N^{*}} - y_{N}|$"]
    
    n_feats = 4
    
    labels_stream_extrs = np.arange(1, n_feats+1)
    
    graphics.plot_log_multi(err_extrs_y[:n_feats], Ns, title_extrs_y, labels_stream_extrs, ax_labels_y)
    graphics.plot_log_multi(err_extrs[:n_feats], Ns, title_extrs, labels_stream_extrs, ax_labels_stream)

    
    # title_left = "Error to $N^{*}=$%d in $y$ of stream-saddle on left boundary"%N_max
    # title_right = "Error to $N^{*}=$%d in $y$ of stream-saddle on right boundary"%N_max
    # ax_labels_saddle = ["N", "$|y_{N^{*}} - y_{N}|$"]
    
    # labels_stream_left = np.arange(1, n_feats+1)
    # labels_stream_right = np.arange(1, n_feats+1)
    
    # graphics.plot_log_multi(err_left[:n_feats], Ns, title_left, labels_stream_left, ax_labels_saddle)
    # graphics.plot_log_multi(err_right[:n_feats], Ns, title_right, labels_stream_right, ax_labels_saddle)

#------------------------------------------------------------------------------

def bfs_compare_N(Ns, N_max):
    
    tru_xs, tru_ys = get_attatchments(N_max)
    
    M = len(Ns)
    a = len(tru_xs)
    b = len(tru_ys)

    err_xs = np.zeros((M, a)) # error in stream extrema along y=0
    err_ys = np.zeros((M, b)) # error in stream extrema along x=x_step
    
    for i in range(M):
        N = Ns[i]
        
        N_xs, N_ys = get_attatchments(N)

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




def plot_bfs_compare_N(Ns, N_max):
    err_xs, err_ys = bfs_compare_N(Ns, N_max)
    print('reattatch', err_xs) 
    print('detatch', err_ys)
    title_xs = "Error to $N^{*}=$%d in reattatchment point"%N_max
    title_ys = "Error to $N^{*}=$%d in detachment point"%N_max
    ax_labels_stream = ["N", "$|\psi_{N^{*}} - \psi_{N}|$"]
    
    n_feats = 4
    
    labels_stream = np.arange(1, n_feats+1)
    
    graphics.plot_log_multi(err_xs[:n_feats], Ns, title_xs, labels_stream, ax_labels_stream)
    graphics.plot_log_multi(err_ys[:n_feats], Ns, title_ys, labels_stream, ax_labels_stream)





    
