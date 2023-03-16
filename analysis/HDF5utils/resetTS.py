#!/usr/bin/env python

###################################################################
# This script works for PSI applications.
# It resets the attribute values to 0 in the concentration group
# corresponding to the last saved time step, including the current 
# time, but except for the surface position, and the delta time 
# set to 1.0e-12 by default.
# The concentrations themselves are not modified. 
# This is useful to use this network file as pre-implanted case,
# starting a new phase of irradiation.
###################################################################

import h5py

## Open the file in which to change the value
f = h5py.File('/home/sophie/Workspace/xolotl-dynamic-build/dynamicGrid/network_dynamic.h5', 'r+')

## Get the last timestep
concGroup = f['concentrationsGroup']
lastTS = concGroup.attrs['lastTimeStep']
lastLoop = concGroup.attrs['lastLoop']

## Open the corresponding group
groupName ='concentrationsGroup/concentration_' + str(lastLoop) + '_' + str(lastTS)
concGroup = f[groupName]

## Get all the attributes
keys = concGroup.attrs.keys()

## Loop on them to reset them
for attrName in keys:
    ## Skip iSurface
    if (attrName == 'iSurface'): continue
    ## deltaTime is set to 1.0e-12 by default
    if (attrName == 'deltaTime'):
         concGroup.attrs[attrName] = 1.0e-12
         continue
    ## reset to 0
    concGroup.attrs[attrName] = 0
