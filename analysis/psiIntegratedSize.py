#!/usr/bin/env python

import numpy as np
import math
import matplotlib.pyplot as plt
from   pylab import *
import h5py
from matplotlib.colors import LogNorm

## Set the Zero
zero = 1.0e-20

## Select the timestep we want to read from
timestep = 9

## Set the maximum size of helium/hydrogen
maxSize = 401

## Open the file (the network/restart HDF5 file)
f = h5py.File('/home/sophie/Workspace/xolotl-plsm-build/script/xolotlStop.h5', 'r')

## Open the concentration group
groupName ='concentrationsGroup/concentration_' + str(timestep)
concGroup = f[groupName]

## Read the concentration and index datasets
concDset = concGroup['concs']
indexDset = concGroup['concs_startingIndices']

## Read the surface position and time at the chosen time step
surfacePos = concGroup.attrs['iSurface']
time = concGroup.attrs['absoluteTime']

## Read the grid to know which grid point is which depth
gridDset = f['headerGroup/grid']

## Read how many normal and super clusters there are
networkGroup = f['networkGroup']
totalSize = networkGroup.attrs['totalSize']

## Create the mesh and data array
heArray = np.empty([maxSize+1])
dArray = np.empty([maxSize+1])
tArray = np.empty([maxSize+1])
hArray = np.empty([maxSize+1])
x = np.empty([maxSize+1])
for i in range(0, maxSize+1):
    heArray[i] = zero
    dArray[i] = zero
    tArray[i] = zero
    hArray[i] = zero
    x[i] = i

## Loop on the grid
for j in range(len(indexDset)-1):
    ## Get the size of the grid point
    dgrid = gridDset[j] - gridDset[j-1]
    ## Loop on the concentrations
    for i in range(indexDset[j], indexDset[j+1]):
        ## Skip the moments for now
        if (int(concDset[i][0]) > totalSize - 1): continue
        ## Get the cluster bounds
        groupName = str(concDset[i][0])
        clusterGroup = networkGroup[groupName]
        bounds = clusterGroup.attrs['bounds']
        if (bounds[8] > 0): continue # I case
        ## Loop on them
        for he in range(bounds[0], bounds[1]+1):
            for d in range(bounds[2], bounds[3]+1):
                for t in range(bounds[4], bounds[5]+1):
                    for v in range(bounds[6], bounds[7]+1):
                        ## Fill the arrays
                        heArray[he] = heArray[he] + concDset[i][1] * dgrid
                        dArray[d] = dArray[d] + concDset[i][1] * dgrid
                        tArray[t] = tArray[t] + concDset[i][1] * dgrid
                        hArray[d+t] = hArray[d+t] + concDset[i][1] * dgrid

## Create plots
fig = plt.figure()
title = 'Concentrations at t = ' + str(time) + ' s'
fig.suptitle(title,fontsize=22)
hePlot = plt.subplot(111)

## Plot the data
hePlot.plot(x, heArray, linewidth=7, color="magenta", label="He")
hePlot.plot(x, dArray, linewidth=6, color="orange", label="D")
hePlot.plot(x, tArray, linewidth=5, color="blue", label="T")
hePlot.plot(x, hArray, linewidth=4, color="cyan", label="H = D+T")
hePlot.set_xlabel("Cluster Size",fontsize=22)
hePlot.set_ylabel("Integrated Concentration (# / nm2)",fontsize=22)
hePlot.set_xlim([1, maxSize])
hePlot.set_ylim([1.0e-10, 10.0])
hePlot.set_yscale('log')
hePlot.tick_params(axis='both', which='major', labelsize=20)

## Plot the legend
l1 = hePlot.legend(loc=9)
setp(l1.get_texts(), fontsize=25)

## Show the plots
plt.show()
