#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 22 14:59:46 2024

@author: sarahdennis
"""
import csv
import time
import numpy as np

from scipy.sparse.linalg import bicgstab
from scipy.sparse.linalg import splu

from scipy.interpolate import interpn
from scipy.signal import argrelextrema as relEx

import graphics
import stokes_solvers as solver

bicgstab_rtol = 1e-8

plot_mod = 50
write_mod = 50
error_mod = 100

import stokes_examples as examples

#------------------------------------------------------------------------------
def run_new(N, iters):
    tri = examples.biswasEx(N)
    # tri = examples.zeroReynEx(N)

    nm = tri.Nx * tri.Ny

    u_init = np.zeros(nm)
    v_init = np.zeros(nm)
    psi_init = np.zeros(nm)
    past_iters = 0

    # u, v, psi = run(tri, u_init, v_init, psi_init, iters, past_iters)
    u, v, psi = run_spLU(tri, u_init, v_init, psi_init, iters, past_iters)
    
    psi = solver.unmirror_boundary(tri, psi)
    write_solution(tri.filename+".csv", nm, u, v, psi, iters)
                                                                                                                                                                                                                                                                             
def run_load(N, iters):

    tri = examples.biswasEx(N)
    # tri = examples.zeroReynEx(N)

    u, v, psi, past_iters = read_solution(tri.filename+".csv", tri.Nx*tri.Ny)
    
    # u, v, psi = run(tri, u, v, psi, iters, past_iters)
    # u, v, psi = run_LU(tri, u, v, psi, iters, past_iters)
    u, v, psi = run_spLU(tri, u, v, psi, iters, past_iters)
    
    psi = solver.unmirror_boundary(tri, psi)
    write_solution(tri.filename+".csv", tri.Nx*tri.Ny, u, v, psi, iters+past_iters)


def load_scale(N_load, N_new):
    # tri_load = triangle(x0, xf, y0, yf, U, Re, N)
    tri_load = examples.biswasEx(N_load) #slope 4 example
    
    
    points_load = (tri_load.ys, tri_load.xs)
    
    u_load, v_load, psi_load, past_iters = read_solution(tri_load.filename+".csv", tri_load.Ny*tri_load.Nx)
    u_load_2D = u_load.reshape((tri_load.Ny,tri_load.Nx), order='F')
    v_load_2D = v_load.reshape((tri_load.Ny,tri_load.Nx), order='F')
    psi_load_2D = psi_load.reshape((tri_load.Ny,tri_load.Nx), order='F')


    tri_scale = examples.biswasEx(N_new)

    points_scale = np.meshgrid(tri_scale.ys, tri_scale.xs)
    
    u_scaled_2D = interpn(points_load, u_load_2D, tuple(points_scale), method='linear')
    v_scaled_2D = interpn(points_load, v_load_2D, tuple(points_scale), method='linear')
    psi_scaled_2D = interpn(points_load, psi_load_2D, tuple(points_scale), method='linear')

    u_scaled = u_scaled_2D.ravel()
    v_scaled = v_scaled_2D.ravel()
    psi_scaled = psi_scaled_2D.ravel()

    write_solution(tri_scale.filename+".csv", tri_scale.Ny*tri_scale.Nx, u_scaled, v_scaled, psi_scaled, 0)
    plot_load(N_new)
    
def plot_load(N):
    # tri = triangle(x0, xf, y0, yf, U, Re, N)
    # tri = biswasEx(N)
    tri = examples.zeroReynEx(N)
    u, v, psi, past_iters = read_solution(tri.filename+".csv", tri.Nx * tri.Ny)

    make_plots(tri, u, v, psi, past_iters)

#------------------------------------------------------------------------------
def read_solution(filename, nm):
    u = np.zeros(nm)
    v = np.zeros(nm)
    psi = np.zeros(nm)
    with open(filename, newline='') as file:
        reader = csv.reader(file)

        for i in range(nm):
            line = next(reader)
            ui, vi, psii = line[0].split(' ')
            u[i] = float(ui)
            v[i] = float(vi)
            psi[i] = float(psii)
        past_iters = int(next(reader)[0])
        file.close()
    return u, v, psi, past_iters

def write_solution(filename, nm, u, v, psi, iters):

    with open(filename, 'w', newline='') as file:
        writer = csv.writer(file, delimiter=' ')
        for i in range(nm):
            writer.writerow([u[i], v[i], psi[i]])
        
        writer.writerow([iters])
        print("  saved csv")
        file.close()    

#------------------------------------------------------------------------------
def run(tri, u, v, past_psi, iters, past_iters):
    
    M = solver.Dpsi_linOp(tri)
    
    for i in range(iters): 
        
        t0 = time.time()
        rhs = solver.update_rhs(tri, u, v)
        
        psi, exit_flag = bicgstab(M, rhs, tol=bicgstab_rtol)
        
        u, v = solver.uv_approx(tri, u, v, psi)
        
        psi = solver.mirror_boundary(tri, psi)

        tf = time.time()
        if i % error_mod == 0: 
            past_psi = solver.unmirror_boundary(tri, past_psi)
            psi = solver.unmirror_boundary(tri, psi)
            
            err_i = np.max(np.abs(psi - past_psi))
            
            print("k=%d of %d"%(i+past_iters+1, iters+past_iters))
            print("  time: %.3f s"%(tf-t0))
            print("  error: %.5e psi"%err_i)
            
            if err_i < bicgstab_rtol:
                print('warning: lower bicgstab_rtol')
            
        if i % plot_mod == 0:
            psi = solver.psi_unmirror_boundary(tri, psi)
            make_plots(tri, u, v, psi, i+1 + past_iters)
        
        if i % write_mod == 0:
            psi = solver.unmirror_boundary(tri, psi)
            write_solution(tri.filename, tri.Nx*tri.Ny, u, v, psi, i+1+past_iters)


        past_psi = psi

    return u, v, psi

def run_spLU(tri, u, v, past_psi, iters, past_iters):
    
    M = solver.Dpsi_cscmatrixBuild(tri)
    LU = splu(M)
    
    for i in range(iters): 
        
        t0 = time.time()
        rhs = solver.update_rhs(tri, u, v)

        psi = LU.solve(rhs)
        
        u, v = solver.uv_approx(tri, u, v, psi)
        
        psi = solver.mirror_boundary(tri, psi)
        
        tf = time.time()
        
        if i % error_mod == 0: 
            past_psi = solver.unmirror_boundary(tri, past_psi)
            psi = solver.unmirror_boundary(tri, psi)
            
            err_i = np.max(np.abs(psi - past_psi))
            
            print("k=%d of %d"%(i+past_iters+1, iters+past_iters))
            print("  time: %.3f s"%(tf-t0))
            print("  error: %.5e psi"%err_i)
            
        if i % plot_mod == 0:
            psi = solver.unmirror_boundary(tri, psi)
            make_plots(tri, u, v, psi, i+1 + past_iters)
        
        if i % write_mod == 0:
            psi = solver.unmirror_boundary(tri, psi)
            write_solution(tri.filename+".csv", tri.Nx*tri.Ny, u, v, psi, i+1+past_iters)


        past_psi = psi

    return u, v, psi



#------------------------------------------------------------------------------
# PLOTTING 
#------------------------------------------------------------------------------
def make_plots(tri, u, v, stream, iters):
    n = tri.Nx
    m = tri.Ny

# Grid domain
    xs = tri.xs
    ys = tri.ys

# Stream: Psi(x,y) heat & contour
    stream_2D = stream.reshape((m,n))
    ax_labels = ['$\psi(x,y)$ : $u = \psi_y$, $v = \psi_x$', '$x$', '$y$']
    title = 'Stream ($N=%d$, $k=%d$)'%(tri.N, iters)
    graphics.plot_contour_heat(stream_2D, xs, ys, title, ax_labels)
      
#  Velocity: (U, V)  streamplot
    u_2D = u.reshape((m,n))
    v_2D = v.reshape((m,n))
    
    ax_labels = ['$|(u,v)|_2$','$x$', '$y$']
    title = 'Velocity ($N=%d$, $k=%d$)'%(tri.N, iters)
    ax_labels = ['$\psi(x,y)$ : $u = \psi_y$, $v = \psi_x$','$x$', '$y$']
    graphics.plot_stream_heat(u_2D, v_2D, xs, ys, stream_2D, title, ax_labels)
    # plot_stream(u_2D, v_2D, xs, ys, title, ax_labels)

#  Vorticity: w = vx - uy heat & contour
    uy_2D = np.gradient(u_2D, tri.dx, axis=0)
    vx_2D = np.gradient(v_2D, tri.dx, axis=1)
    w = np.zeros((m,n))
    for j in range(m):
        for i in range(n):   
            w[j,i] = vx_2D[j,i] - uy_2D[j,i]
    ax_labels = ['$\omega(x,y) = -( \psi_{xx} + \psi_{yy})$', '$x$', '$y$']
    title = 'Vorticity ($N=%d$, $k=%d$)'%(tri.N, iters)
    graphics.plot_contour_heat(w, xs, ys, title, ax_labels)


#------------------------------------------------------------------------------
# Critiacal points 
#------------------------------------------------------------------------------

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

def write_criticals(N):
    tri = examples.biswasEx(N)
    u, v, psi, past_iters = read_solution(tri.filename+".csv", tri.Nx * tri.Ny)
    
    left, right = get_boundary(tri, psi)
    
    crits_filename = 'crits_' + tri.filename +".csv"
    with open(crits_filename, 'w', newline='') as file:
        writer = csv.writer(file, delimiter=' ')
    
        # print("Psi sign changes...")
        writer.writerow('xy')
        sign_ref = 0
        for (x, y, p) in left:
            sign_new = np.sign(p)
            if sign_new != 0 and sign_new != sign_ref:
                sign_ref = sign_new
                # print("(x:%.5f, y:%.6f)"%(x,y))
                writer.writerow([x,y])
                
        writer.writerow('xy')
        sign_ref = 0
        for (x, y, p) in right:
            sign_new = np.sign(p)
            if sign_new != 0 and sign_new != sign_ref:
                sign_ref = sign_new
                # print("(x:%.5f, y:%.6f)"%(x,y))
                writer.writerow([x,y])

    
        # print("Psi center-line extrema...")
        writer.writerow('xyp')
        center = get_center(tri, psi)
        max_inds = relEx(center[:,2], np.greater)[0]
        min_inds = relEx(center[:,2], np.less)[0]
        
    
        for i in max_inds:
            x,y,p = center[i]
            # print("(x:%.1f, y:%.6f) p=%.5e"% (x,y,p))
            writer.writerow([x,y,p])
        
        writer.writerow('xyp')
        for i in min_inds:
            x,y,p = center[i]
            # print("(x:%.1f, y:%.6f) p=%.5e"% (x,y,p))
            writer.writerow([x,y,p])

    make_plots(tri, u, v, psi, past_iters)
