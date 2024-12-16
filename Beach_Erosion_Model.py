#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 12 10:21:33 2024

@author: westanlubin
"""

import numpy as np
import matplotlib.pyplot as plt
import sys

### -------- INITIAL CONDITIONS -------------

# Load elevation data from an ASCII file (e.g., a DEM). The first 6 rows contain metadata.
ascii_grid = np.loadtxt("alabeach.asc", dtype='float', skiprows=6)

# Replace missing data values (-9999) with NaN to handle them correctly in calculations.
ascii_grid[ascii_grid == -9999] = np.nan

# Cap elevation values above 2 at 2, representing some threshold for the dataset.
ascii_grid[ascii_grid >= 2] = 2

# Load the ASCII headers (metadata), which contain grid size, resolution, and coordinates.
ascii_headers = np.loadtxt("alabeach.asc", max_rows=6, dtype='str')
n_long = ascii_headers[0, 1].astype(int)  # Number of columns (longitude)
n_lat = ascii_headers[1, 1].astype(int)   # Number of rows (latitude)
dxy = ascii_headers[4, 1].astype(float)   # Grid cell size (resolution)
xllcorner = ascii_headers[2, 1].astype(float)  # X-coordinate of the lower-left corner
yllcorner = ascii_headers[3, 1].astype(float)  # Y-coordinate of the lower-left corner

# Remove a 2-cell boundary around the grid (likely for numerical stability or cropping).
ascii_grid = ascii_grid[2:-2, 2:-2]

# Update grid dimensions after trimming.
n_lat, n_long = ascii_grid.shape

### GRID SETUP
nx = n_lat  # Number of rows (latitude)
ny = n_long  # Number of columns (longitude)

dx = 1  # Grid spacing in x (meters)
dy = 1  # Grid spacing in y (meters)

# Create coordinate arrays for the grid.
x = np.arange(0, dxy * nx, dxy) + xllcorner  # X-coordinates of grid cells
y = np.arange(0, dxy * ny, dxy) + yllcorner  # Y-coordinates of grid cells
X, Y = np.meshgrid(x, y, indexing='ij')  # Create a 2D grid for plotting.

# Total number of grid nodes.
nodes = n_long * n_lat

### TIME STEP PARAMETERS
dt = 5  # Time step in years

# Initialize uniform velocity across the grid.
u = np.full((nx, ny), 0.01)  # Velocity (m²/year).

# Modify the velocity on the rightmost 20% of the grid to simulate shoreline hardening effects.
right_strip = int(0.8 * nx)  # Adjust to cover the last 20% of columns

# Assign higher velocity values in this strip
for i in range(right_strip, ny):  # Loop through the columns in the strip
    u[i] = 0.04

# Flatten the velocity grid to match the elevation data's 1D format.
u_flat = u.flatten()

### INITIAL ELEVATION SETUP
# Flatten the ASCII grid to create a 1D elevation array for calculations.
z = ascii_grid.flatten()

### STABILITY CHECKS
# Courant numbers in x and y directions to ensure numerical stability.
cx = dt * np.max(u) / dx
cy = dt * np.max(u) / dy

# Check for invalid values (e.g., NaN in the elevation grid).
if np.sum(np.isnan(ascii_grid)) > 0:
    print('NaN values detected in the elevation grid.')
    sys.exit()

# Ensure stability conditions are met; exit if unstable.
if cx > 1:
    print('Stability condition violated: cx > 1')
    sys.exit()
if cy > 1:
    print('Stability condition violated: cy > 1')
    sys.exit()

### MATRIX A FORMATION FOR TIME EVOLUTION
# Create the A matrix, which determines how elevation evolves over time.
A = np.zeros((nx * ny, nx * ny))

for i in range(nx):  # Loop through rows.
    for k in range(ny):  # Loop through columns.
        ik = i * ny + k  # Flattened index of the current node.
        ## BOUNDARY CONDITIONS
        if i == 0 or k == 0:  # For the first row/column, no change.
            A[ik, ik] = 1
        else:
            ## MATRIX COEFFICIENTS
            A[ik, ik] = 1 - cx - cy  # Contribution from the current cell.
            A[ik, (i - 1) * ny + k] = cx  # Contribution from the cell above.
            A[ik, i * ny + k - 1] = cy  # Contribution from the cell to the left.

### PLOT INITIAL CONDITIONS
# Reshape elevation data for plotting.
Z = z.reshape(X.shape)

# Plot initial elevation conditions using a color map.
fig, ax = plt.subplots(1, 1)
c1 = ax.pcolormesh(X, Y, Z, cmap='viridis')
fig.colorbar(c1)
ax.set_xlabel('Distance (m)')
ax.set_ylabel('Depth (m)')
ax.set_title('Initial Elevation Conditions')

# Create a 3D surface plot of the initial elevation.
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot_surface(X, Y, Z)
ax.set_xlabel('Distance (m)')
ax.set_ylabel('Depth (m)')
ax.set_title('Initial Elevation Surface')

### RUN TIME LOOP
# Simulate elevation change over time.
totaltime = 100  # Total simulation time (years).
time = 0
while time <= totaltime:
    newz = np.dot(A, z)  # Compute the new elevation using matrix multiplication.
    z[:] = newz  # Update the elevation.
    time += dt  # Increment the time.

### PLOT FINAL CONDITIONS
# Reshape the final elevation data for plotting.
Z = z.reshape(X.shape)

# Create a 3D surface plot of the final elevation.
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot_surface(X, Y, Z)
ax.set_xlabel('Distance (m)')
ax.set_ylabel('Depth (m)')
ax.set_title('Final Elevation Surface')

# Create a 2D color map of the final elevation.
fig2, ax2 = plt.subplots(1, 1)
c = ax2.pcolormesh(X, Y, Z)
ax2.set_title('Final Elevation Conditions')
colorbar = fig2.colorbar(c, ax=ax2)
colorbar.set_label('Elevation (m)')

### PLOT SPATIALLY VARYING VELOCITIES
# Reshape velocity data for plotting.
U = u.reshape(X.shape)

# Create a 2D color map of the spatially varying velocity.
fig, ax = plt.subplots(1, 1)
c = ax.pcolormesh(X, Y, U, cmap='plasma')  # Use a colormap to show velocity variations.
fig.colorbar(c, ax=ax, label='Velocity (m²/year)')
ax.set_title('Spatially Varying Velocity')
ax.set_xlabel('X Coordinate (m)')
ax.set_ylabel('Y Coordinate (m)')

plt.show()

### SAVE ASCII OUTPUT
# Save the final elevation data to an ASCII file for GIS visualization.
header = f'NCOLS {n_long} \nNROWS {n_lat} \nxllcorner {xllcorner} \nyllcorner {yllcorner} \ncellsize {dxy} \nNODATA_value -9999'
np.savetxt('InitialBeach_elev.asc', z, header=header, comments='')






