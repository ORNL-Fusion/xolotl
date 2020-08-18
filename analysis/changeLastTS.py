#!/usr/bin/env python

import h5py

## Set the new TS value to use
newTS = 7308

## Open the file in which to change the value
f = h5py.File('/home/sophie/Workspace/xolotl-dynamic-build/dynamicGrid/network_dynamic.h5', 'r+')
## Modify the last time step saved in the file
concGroup = f['concentrationsGroup']
concGroup.attrs['lastTimeStep'] = newTS
