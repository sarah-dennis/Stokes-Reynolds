# -*- coding: utf-8 -*-
"""
Created on Fri May 24 08:38:16 2024

@author: sarah
"""
import csv
import numpy as np
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
        err = int(next(reader)[0])
        file.close()
    return u, v, psi, err

def write_solution(ex, u, v, psi, err):
    nm = ex.Nx * ex.Ny
    filename = ex.filestr + ".csv"
    with open(filename, 'w', newline='') as file:
        writer = csv.writer(file, delimiter=' ')
        for i in range(nm):
            writer.writerow([u[i], v[i], psi[i]])
        
        writer.writerow([err])
        file.close()    
