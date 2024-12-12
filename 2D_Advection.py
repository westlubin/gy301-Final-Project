#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 12 10:47:31 2024

@author: westanlubin
"""
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

### INITIAL CONDITIONS
nx = 15
ny = 20

dx = 2 # meters
dy = 2 # meters

x = np.arange(0,nx*dx,dx) # created 1-D array of x positions
y = np.arange(0,ny*dy,dy) # created 1-D array of x positions
X,Y = np.meshgrid(x,y,indexing = 'ij') # created 2-D coordinate system for plotting

""" 
capital letters = 2d array
lowercase = 1d array 
"""

u = 0.02 # m^2/yr

dt = 5 # years

# example will use topography
Z = np.random.random((nx,ny))*100
z=Z.flatten() # gives the 1-D flattened array, flattened by rows

### STABILITY CHECK
cx = dt*u/dx
cy = dt*u/dy

import sys
if cx>1:
    print('x is unstable')
    sys.exit()

if cy>1:
    print('y is unstable')
    sys.exit()

### CREATING THE A MATRIX

A = np.zeros((nx*ny, nx*ny))

for i in range(nx):
    for k in range(ny):
        ik = i*ny + k
        ## -----BOUNDRY CONDITIONS----- ##
        if i==0:
            A[ik,ik] = 1 # no change
        elif k==0:
            A[ik,ik] = 1
        else:
            ## -----MATRIX COEFFICIENT----- ##
            A[ik,ik] = 1-cx - cy
            A[ik,(i-1)*ny + k] = cx
            A[ik,i*ny + k - 1] = cy
#print(A)

# method 1 - use a surface plot
fig = plt.figure()
ax = fig.add_subplot(111, projection = '3d')
ax.plot_surface(X,Y,Z)

#method 2 - uses pcolormesh
fig2, ax2 = plt.subplots(1,1)
c = ax2.pcolormesh(X,Y,Z)
ax2.set_title('initial conditions')
colorbar = fig2.colorbar(c, ax=ax2)
colorbar.set_label('Elevation')
# Add legend or color bar

### RUNNING TIME
totaltime = 1000
time = 0
while time <= totaltime: 
    newz = np.dot(A,z)
    z[:] = newz
    time += dt

### PLOTTING THE FINAL TOPOGRAPHY
Z = z.reshape(X.shape)
fig = plt.figure()
ax = fig.add_subplot(111,projection = '3d')
ax.plot_surface(X,Y,Z)

fig2, ax2 = plt.subplots(1,1)
c = ax2.pcolormesh(X,Y,Z)
ax2.set_title('final conditions')
colorbar = fig2.colorbar(c, ax=ax2)
colorbar.set_label('Elevation')