readme.txt for Xolotl results
20160527 Sophie Blondel

The grainBoundariesStudy folder contains parameter files in order 
to reproduce results for the study of small helium clusters behavior near 
the GB of tungsten material.

The Xolotl version used here is revision 1012 of the irregular_grid_branch.

The network file networkInit.h5 contains a reduced network up to V = 250, 
corresponding to 1581 clusters. The grid is a 2D grid with 100x50 grid 
points, the step size is variable in the x direction (depth) and 1.0 nm
in the y direction. The void portion is 41.0 %.

The parameter file param.txt includes two Y GB at 12.0 and 38.0 nm.

To obtain the concentration of helium as a function of depth you have to use 
the -helium_conc PETSc option and keep the values either:
- on the GBs: +/- 2nm around 12.0 and 38.0
- or off of the GBs: below 2.0, from 23.0 to 27.0, and above 47.0 nm

To reproduce the result with mobile clusters on the GBs, you have to 
modify the code manually. In the diffusion handler you should force 
all the values to be 1.0 on the grid.