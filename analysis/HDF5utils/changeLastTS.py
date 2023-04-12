#!/usr/bin/env python

###################################################################
# This script works for any applications.
# It lets the user define a new last time step (line 14) to read
# from when this network file is loaded by Xolotl.
# The new TS number should correspond to an existing concentration
# group.
###################################################################

import h5py

## Set the new TS value to use
newTS = 7308
newLoop = 7308

## Open the file in which to change the value
f = h5py.File('/home/sophie/Workspace/xolotl-dynamic-build/dynamicGrid/network_dynamic.h5', 'r+')
## Modify the last time step saved in the file
concGroup = f['concentrationsGroup']
concGroup.attrs['lastTimeStep'] = newTS
concGroup.attrs['lastLoop'] = newLoop
