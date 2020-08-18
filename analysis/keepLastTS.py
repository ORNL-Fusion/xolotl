#!/usr/bin/env python

import h5py

## Open the file we want to copy from
f = h5py.File('/home/sophie/Data/Xolotl/network/networkTest.h5', 'r')

## Get the last time step saved in the file
concGroup = f['concentrationsGroup']
timestep = concGroup.attrs['lastTimeStep']

## Create the file to copy to
fNew = h5py.File('/home/sophie/Data/Xolotl/network/networkNew.h5', 'a')

## Create the concentration group
concGroupNew = fNew.create_group('concentrationsGroup')

## Set the last time step
concGroupNew.attrs['lastTimeStep'] = timestep

## Copy the last timestep group
groupName ='concentration_' + str(timestep)
concGroup.copy(groupName, concGroupNew)

## Copy the other groups
f.copy('headerGroup', fNew)
f.copy('networkGroup', fNew)
