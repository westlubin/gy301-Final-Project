MODEL PURPOSE

The purpose of this model is to simulate the evolution of an elevation grid representing a 100 x 100-meter area of beach located in the Ala Moana Beach Park of Honolulu, HI.
The model runs for 100 years with spatially varying velocity to represent increased erosion effects due to shoreline hardening.
The model uses finite difference methods to calculate how elevation changes under specific physical conditions, ensuring numerical stability with Courant criteria.

MODEL INPUTS

A digital elevation model in ASCII format (alabeach.asc) representing the initial elevation grid. 
Velocity field (u) specifies the rate of elevation change in space and time. 
Uniform velocity is initialized, and higher velocities are applied to the rightmost 20% of the grid, simulating the effects of shoreline hardening.
Spatial and temporal parameters: grid dimensions (nx,ny) derived from the ASCII file, spatial step for x and y directions (dx,dy), time step 
to determine the temporal resolution of the simulation (dt), and the total period for which the model runs (totaltime). 

MODEL OUTPUTS

2D and 3D plots for the initial and final conditions show how the surface evolved over 100 years. A 2D colormap of the spatially
Varrying velocity field is used for visualization of where the intended impact of shoreline hardening is on the model. 
An updated ASCII file containing the final simulated elevation field after the 100-year time period (FinalBeach_elev.asc).
