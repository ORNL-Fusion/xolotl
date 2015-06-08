readme.txt for Xolotl results
20150608 Sophie Blondel

The modifiedTrapMutationStudy folder contains parameter files in order 
to reproduce results for the study of small helium clusters behavior near 
the surface of tungsten material.

The Xolotl version used here is revision 882 of the irregular_grid_branch.

Two network files are included: 
- network.h5: a full one up to V = 50 (corresponding to 5448 clusters)
- networkReduced.h5: a reduced one up to V = 50 too. For the reduced one only 
bubbles close to He/V = 4.0 are kept above V = 12, so that the total number 
of clusters is only 581.
They both contain a one dimensional grid of 59 point in order to reach a depth 
of 25.5 nm when the option "regularGrid=no" is used in the parameter file.

The parameter files use network.h5, a W100 surface, all the processes 
(reaction, diffusion, advection, modified trap-mutation, desorption), and the 
"-helium_retention -helium_cumul" options for outputs are set by default.

Three parameter files are included, setup to run up to a fluence of 10^20 He/m^2 
at 933 K:
- param_4e25.txt: for a flux of 4.0x10^25 He/m^2/s
- param_1.6e26.txt: for a flux of 1.6x10^26 He/m^2/s
- param_5e27.txt: for a flux of 5.0x10^27 He/m^2/s
