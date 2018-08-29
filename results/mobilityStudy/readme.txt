readme.txt for Xolotl results
20170514 Sophie Blondel

The mobilityStudy folder contains parameter files in order 
to reproduce results for the study of small helium-vacancy clusters mobility 
in tungsten material.

The Xolotl version used here is revision 1171 of the trunk.

One network file in included: 
- networkMobile.h5: a reduced network up to V = 50. For the reduced one only 
bubbles close to He/V = 4.0 are kept above V = 12, so that the total number 
of clusters is only 1552, this cut follow the phase-cut method. In this file
a few HeV clusters have their diffusion factors set to some values.
It contains a one dimensional grid of 200 points in order to reach a 
depth of 49.5 nm when the option "regularGrid=no" is used in the parameter 
file associated to a void portion of 58.5%.

The parameter file uses networkMobile.h5, a W100 surface, all the processes 
(reaction, diffusion, advection, modified trap-mutation, desorption, moving 
surface), and the "-helium_retention -helium_conc" options for outputs 
are set by default. It is set-up to run up to a fluence of 10^22 He/m^2 
at 1000 K.
