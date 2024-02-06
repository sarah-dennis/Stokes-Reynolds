#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 22 14:59:46 2024

@author: sarahdennis
"""
import numpy as np

from scipy.sparse.linalg import LinearOperator
from scipy.sparse.linalg import bicgstab

import graphics as graph

x0 = 0
xL = 1

y0 = 0
yL = 2

slope = yL/(xL/2)
U = 1

Re = 1

N = 100
h = 1/N

n = xL*N + 1
m = yL*N + 1


# concatenated wedge slider bearings 

# write to the new Biswas email

# if psi is 0, and boundary velocity's given, for (14) and (15)


#rhs = c0*A - c1*(B - C)    
c0 = 3*h
c1 = 0.5 * h**2 * Re
A, B, C = np.zeros(3)
# A = u_S - u_N + v_E - v_W
# B = v_C * (u_E + u_W + u_N + u_S)
# C = u_C * (v_E + v_W + v_N + v_S)

def make_rhs(u, v): #
    rhs = np.zeros(m*n)

    for k in range(n*m):
        i = k % n
        j = k//n
        
        x = h*i
        y = h*j
       
        # boundary
        if i == 0 or j == 0 or i == n-1 or j == m-1:
            rhs[k] = 0
            # TODO even for j=0 where ux=U ?
          
        # dead zones 
        elif y <= -slope * x + yL or y <= slope * x - yL:
            rhs[k] = 0
        
        # interior
        else:
            # (u,v) at 9 point stencil
            #k = j*n + i
            u_C = u[k]
            v_C = v[k]
            
            # North (i, j+1)
            if j+1 == m-1: 
                u_N = U
                v_N = 0
            else:
                u_N = u[(j+1)*n + i]
                v_N = v[(j+1)*n + i]
                
            # East (i+1, j)
            x_E = h*(i+1)
            if i+1 == n-1 or y <= slope * x_E - yL: 
                u_E = 0
                v_E = 0
            else:
                v_E = v[j*n + i+1]
                u_E = u[j*n + i+1]
            
            # South (i, j-1)
            y_S = h*(j-1)
            if j-1 == 0 or y_S <= slope*x - yL or y_S <= -slope*x + yL: 
                v_S = 0
                u_S = 0
            else:
                v_S = v[(j-1)*n + i]
                u_S = u[(j-1)*n + i]
                
            # West (i-1, j)  
            x_W = (i-1)*h
            if i-1 == 0 or y <= -slope * x_W - yL: 
                v_W = 0
                u_W = 0
            else:
                v_W = v[j*n + (i-1)]
                u_W = v[j*n + (i-1)]
                
            A = u_S - u_N + v_E - v_W
            B = v_C * (u_E + u_W + u_N + u_S)
            C = u_C * (v_E + v_W + v_N + v_S)

            rhs[k] = c0 * A + c1 * (B - C)
            
    return rhs

#
# new_u[k] = c2 * (psi_N - psi_S) - c3 * (u_N + u_S)
# new_v[k] = -c2 * (psi_E - psi_W) - c3 * (v_E + v_W)
c2 = 3/(4*h)
c3 = 1/4

def uv_approx(u, v, psi):
    new_u = np.zeros(n*m)
    new_v = np.zeros(n*m)

    for k in range(n*m):
        i = k % n
        j = k//n
        
        x = h*i
        y = h*j
    
        # y=yL boundary
        if j == m-1: 
            new_u[k] = U
            new_v[k] = 0 
            
        # other side boundaries
        elif i == 0 or j == 0 or i == n-1:
            new_u[k] = 0
            new_v[k] = 0 
            
        # dead zones
        elif y <= -slope * x + yL or y <= slope * x - yL:
            
            new_u[k] = 0
            new_v[k] = 0 
                
        else: #interior
            
            # (u,v) at 4 point stencil

            # North (i, j+1)
            if j+1 == m-1: 
                u_N = U
                psi_N = 0
            else:
                u_N = u[(j+1)*n + i]
                psi_N = psi[(j+1)*n + i]
                
            # East (i+1, j)
            x_E = h*(i+1)
            if i+1 == n-1 or y <= slope * x_E - yL: 
                v_E = 0
                psi_E = 0
            else:
                v_E = v[j*n + i+1]
                psi_E = psi[j*n + i]
            
            # South (i, j-1)
            y_S = h*(j-1)
            if j-1 == 0 or y_S <= slope*x - yL or y_S <= -slope*x + yL: 
                u_S = 0
                psi_S = 0
            else:
                u_S = u[(j-1)*n + i]
                psi_S = psi[(j-1)*n + i]
                
            # West (i-1, j)  
            x_W = h*(i-1)
            if i-1 == 0 or y <= -slope * x_W - yL: 
                v_W = 0
                psi_W = 0
            else:
                v_W = v[j*n + (i-1)]
                psi_W = psi[j*n + i-1]
   
            new_u[k] = c2 * (psi_N - psi_S) - c3 * (u_N + u_S)
            new_v[k] = -c2 * (psi_E - psi_W) - c3 * (v_E + v_W)
    
    return new_u, new_v

class linOp(LinearOperator):

    def __init__(self):
        #  global vars: Nx, Ny, h, slope, yL
        
        self.shape = ((n*m, n*m))        
        self.dtype = np.dtype('f8')
        self.mv = np.zeros(n*m) 
        
    def _matvec(self, v): # v:= psi[i*n + j]  

        for k in range(n*m):
            i = k%n
            j = k//n
            
            x = h*i
            y = h*j


            # outer boundary
            if i == 0 or i == n-1 or j == 0 or j == m-1:
                self.mv[k] = 0
    
            # dead zones: 
            elif y <= -slope * x + yL or y <= slope * x - yL:
                self.mv[k] = 0

            #interior   
            else: 
                            
                v_C = v[k]
        
                # psi[k] at 9 point stencil
                
                #North (i, j+1)
                y_N = h*(j+1)
                if j+1 == m-1: 
                    v_N = 0
                else:
                    v_N = v[(j+1)*n + i]
                    
                #East (i+1, j)
                x_E = h*(i+1)
                if i+1 == n-1 or y <= slope * x_E - yL:
                    v_E = 0
                else:
                    v_E = v[j*n + (i+1)]
                        
                #South (i, j-1)
                y_S = h*(j-1)
                if j-1 == 0 or y_S <= slope*x - yL or y_S <= -slope*x + yL: 
                    v_S = 0
                else:
                    v_S = v[(j-1)*n + i]
                
                #West (i-1,j)
                x_W = h*(i-1)
                if i-1 == 0 or y <= -slope * x_W - yL: 
                    v_W = 0
                else:
                    v_W = v[j*n + (i-1)]
                    
                #NorthEast (i+1, j+1)
                if i+1 == n-1 or j+1 == m-1 or y_N <= slope * x_E - yL: 
                    v_NE = 0
                else:
                    v_NE  = v[(j+1)*n + (i+1)]
                
                #NorthWest (i-1, j+1)
                if i-1 == 0 or j+1 == m-1 or y_N <= -slope * x_W - yL:
                    v_NW = 0
                else:
                    v_NW = v[(j+1)*n + (i-1)]
                
                #SouthEast (i+1, j-1)
                if i+1 == n-1 or j-1 == 0 or y_S <= slope*x_E - yL or y_S <= -slope*x_E + yL: 
                    v_SE = 0
                else:
                    v_SE = v[(j-1)*n + i+1]
                
                #SouthWest (i-1, j-1)
                if i-1 == 0 or j-1 == 0 or y_S <= slope*x_W - yL or y_S <= -slope*x_W + yL:
                    v_SW = 0
                else:
                    v_SW = v[(j-1)*n + i-1]

        
                self.mv[k] = v_SW - 8*v_S + v_SE + -8*v_W + 28*v_C -8*v_E + v_NW -8*v_N + v_NE

        return self.mv


M = linOp()

u = np.ones(n*m)
v = np.ones(n*m)

stab_tol = 1e-5
trials = 10

for i in range(trials): 
    print("trial %d of %d"%(i+1, trials))
    rhs = make_rhs(u, v)
    stream, exit_flag = bicgstab(M, rhs, atol=stab_tol)
    u, v = uv_approx(u, v, stream)

#graphing...
xs = np.linspace(x0, xL, n)
ys = np.linspace(y0, yL, m)
u = u.reshape((m,n))
v = v.reshape((m,n))
ax_labels = ['x', 'y']
title = 'velocity stream-plot ($N=%d$, $trials=%d$)'%(N, trials)
graph.plot_stream(u, v, xs, ys, title, ax_labels)











