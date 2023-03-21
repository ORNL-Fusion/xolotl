#!/usr/bin/env python

###################################################################
# This script works for any applications.
# It creates a new file (line 21) containing only the last 
# concentration group and everything else needed to restart Xolotl.
# This is useful when restarting the simulation often and when 
# the other concentration groups won't be needed for data analysis.
###################################################################

import h5py

## Open the file we want to copy from
f = h5py.File('/home/sophie/Workspace/xolotl-dynamicGrid-build/script/xolotlStop.h5', 'r')

## Get the last time step saved in the file
concGroup = f['concentrationsGroup']
timestep = concGroup.attrs['lastTimeStep']
lastLoop = concGroup.attrs['lastLoop']

## Create the file to copy to
fNew = h5py.File('/home/sophie/Workspace/xolotl-dynamicGrid-build/script/networkNew.h5', 'a')

## Create the concentration group
concGroupNew = fNew.create_group('concentrationsGroup')

## Set the last time step
concGroupNew.attrs['lastTimeStep'] = timestep
concGroupNew.attrs['lastLoop'] = lastLoop

## Copy the last timestep group
groupName ='concentration_' + str(lastLoop) + '_' + str(timestep)
concGroup.copy(groupName, concGroupNew)

## Copy the other groups
f.copy('networkGroup', fNew)
