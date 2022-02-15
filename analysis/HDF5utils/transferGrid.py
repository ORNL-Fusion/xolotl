#!/usr/bin/env python

###################################################################
# This script works for PSI applications.
# It lets the user increase or decrease the size of the grid by
# a number of grid points defined in line 19  (nPoints can be 
# negative to decrease the size of the grid).
# The script will print two values to be updated in the parameter
# file: new void portion and new number of grid points.
###################################################################

import numpy as np
import math
import matplotlib.pyplot as plt
from   pylab import *
import h5py

## The number of grid points to be added
nPoints = 25

## Open the file
f0 = h5py.File('/home/sophie/Workspace/xolotl-release-build/transferGrid/networkRef.h5', 'r+')
## Read the grid to know which grid point is which depth
gridDset0 = f0['headerGroup/grid']

## Look for the original grid position
for i in range(1,len(gridDset0)):
    if (gridDset0[i] - gridDset0[i-1] > 0.11): break
initSurf = i - 25
## Compute the current void portion
voidPortion = 100.0 * float(initSurf) / float(len(gridDset0))
print('Old void portion:', voidPortion, 'with', len(gridDset0), 'grid points.')
## Compute the new one
voidPortion = 100.0 * float(initSurf + nPoints) / float(len(gridDset0) + nPoints)
print('New void portion:', voidPortion, 'with', len(gridDset0) + nPoints, 'grid points.')

## Add points to the new grid if nPoints>0
newGrid = []
if (nPoints > 0):
    for i in range(nPoints):
        newGrid.append(i * 0.1)
    for i in range(len(gridDset0)):
        newGrid.append(gridDset0[i] + (nPoints) * 0.1)
## Remove them
else:
    for i in range(-nPoints, len(gridDset0)):
        newGrid.append(gridDset0[i] + (nPoints) * 0.1)

## Replace the grid dataset
gridArray = np.array(newGrid, dtype=np.float)
del f0['headerGroup/grid']
dataset = f0.create_dataset('headerGroup/grid', (len(gridArray),), dtype=np.float)
dataset[...] = gridArray

## Open the header group
headerGroup = f0['headerGroup']
## Replace the nx attribute
nx = headerGroup.attrs['nx']
headerGroup.attrs['nx'] = nx + nPoints

## Get the last time step saved in the file
concGroup0 = f0['concentrationsGroup']
timestep = concGroup0.attrs['lastTimeStep']
## Open the concentration group at this time step
groupName ='concentration_' + str(timestep)
subConcGroup0 = concGroup0[groupName]
## Read the concentration and index datasets
concDset0 = subConcGroup0['concs']
indexDset0 = subConcGroup0['concs_startingIndices']
## Replace the surface position
surfacePos0 = subConcGroup0.attrs['iSurface']
subConcGroup0.attrs['iSurface'] = surfacePos0 + nPoints

## Read the last data in the conc dataset, it corresponds to the temperature and its index
tempData = concDset0[len(concDset0)-1]
## Create the new concentrations
newConc = []
if (nPoints > 0):
    for i in range(nPoints):
        tempConc = (tempData[0], tempData[1])
        newConc.append(tempConc)
    for i in range(len(concDset0)):
        tempConc = (concDset0[i][0], concDset0[i][1])
        newConc.append(tempConc)
else:
    for i in range(-nPoints, len(concDset0)):
        tempConc = (concDset0[i][0], concDset0[i][1])
        newConc.append(tempConc)

## Create a specific data type for concentrations
dType = np.dtype([('ConcType.first', np.int32),
                  ('ConcType.second', np.float)])
concArray = np.array(newConc, dtype=dType)
## Replace the concentration dataset
del subConcGroup0['concs']
dataset = subConcGroup0.create_dataset('concs', (len(concArray),), dtype=dType)
dataset[...] = concArray
## Create the new conc indices
newIndex = []
if (nPoints > 0):
    for i in range(nPoints):
        newIndex.append(i)
    for i in range(len(indexDset0)):
        newIndex.append(indexDset0[i] + nPoints)
else:
    for i in range(-nPoints, len(indexDset0)):
        newIndex.append(indexDset0[i] + nPoints)
indexArray = np.array(newIndex, dtype=np.uint32)
del subConcGroup0['concs_startingIndices']
dataset = subConcGroup0.create_dataset('concs_startingIndices', (len(indexArray),), dtype=np.uint32)
dataset[...] = indexArray
