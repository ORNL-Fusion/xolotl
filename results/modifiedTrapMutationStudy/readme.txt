readme.txt for Xolotl results
20150623 Sophie Blondel

The modifiedTrapMutationStudy folder contains parameter files in order 
to reproduce results for the study of small helium clusters behavior near 
the surface of tungsten material.

The Xolotl version used here is revision 894 of the irregular_grid_branch.

Two network files are included: 
- network50.h5: a reduced network up to V = 50. For the reduced one only 
bubbles close to He/V = 4.0 are kept above V = 12, so that the total number 
of clusters is only 581, this cut follow the phase-cut method.
- network250.h5: a reduced network up to V = 250, corresponding to 1581 
clusters.
They both contain a one dimensional grid of 200 points in order to reach a 
depth of 25.5 nm when the option "regularGrid=no" is used in the parameter 
file associated to a void portion of 70.5%.

The parameter file uses network250.h5, a W100 surface, all the processes 
(reaction, diffusion, advection, modified trap-mutation, desorption, moving 
surface), and the "-helium_retention -helium_cumul" options for outputs 
are set by default. It is set-up to run up to a fluence of 10^20 He/m^2 
at 933 K.

Two types of results are included: helium retention ("retention"), as well 
as cumulative helium distribution ("heliumCumul").

The first column of the retention file represents the fluence in He/nm^2 
and the second one is the percentage of retention. In the heliumCumul file 
the first column is the absolute depth in nm (you need to substract 14.1nm
to get the depth relative to the original surface position), and the 
second one is the cumulative quantity of helium at this depth (not 
normalized).

The naming convention for these files is the following: first is the surface 
orientation (W001, W011, W111, and W211), second is the list of processes 
that were used (1: diffusion, 2: advection, 3: modified trap-mutation, 4:
moving surface).

For the cumulative distribution, the information is kept at a fluence of 
~3.3e19 He/m^2, except if "low" is appended to the file name (in that case
it is saved at a fluence of ~3-5e18 He/m^2).
