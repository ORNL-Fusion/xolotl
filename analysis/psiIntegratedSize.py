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
timestep = 0

## Set the maximum size of helium/hydrogen
maxSize = 401

## Open the file (the network/restart HDF5 file)
f = h5py.File('/home/sophie/Data/Xolotl/network/He2DT/network_ITER_BPO_37_2020.h5', 'r')

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

## Read the composition index to know which cluster is what
compDset = f['headerGroup/composition']

## Read how many normal and super clusters there are
networkGroup = f['networkGroup']
normalSize = networkGroup.attrs['normalSize']
superSize = networkGroup.attrs['superSize']

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
        if (int(concDset[i][0]) > len(compDset) - 1): continue
        ## Normal clusters
        if (int(concDset[i][0]) < normalSize):
            ## Get the helium and hydrogen sizes of this cluster
            heSize = compDset[int(concDset[i][0])][0]
            dSize = compDset[int(concDset[i][0])][1]
            tSize = compDset[int(concDset[i][0])][2]
            ## Fill the arrays
            heArray[heSize] = heArray[heSize] + concDset[i][1] * dgrid
            dArray[dSize] = dArray[dSize] + concDset[i][1] * dgrid
            tArray[tSize] = tArray[tSize] + concDset[i][1] * dgrid
            hArray[dSize+tSize] = hArray[dSize+tSize] + concDset[i][1] * dgrid
        ## Super clusters
        else:
            ## Loop on the number of clusters it contains
            groupName = str(concDset[i][0]) + '/heVList'
            heVList = networkGroup[groupName]
            for k in range(0, len(heVList)):
                ## Get the helium and hydrogen sizes of this cluster
                heSize = heVList[k][0]
                dSize = heVList[k][1]
                tSize = heVList[k][2]
                ## Fill the arrays
                heArray[heSize] = heArray[heSize] + concDset[i][1] * dgrid
                dArray[dSize] = dArray[dSize] + concDset[i][1] * dgrid
                tArray[tSize] = tArray[tSize] + concDset[i][1] * dgrid
                hArray[dSize+tSize] = hArray[dSize+tSize] + concDset[i][1] * dgrid

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
