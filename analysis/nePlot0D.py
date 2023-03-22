#!/usr/bin/env python

import numpy as np
import math
import matplotlib.pyplot as plt
from   pylab import *
import h5py
from matplotlib.colors import LogNorm

## Set the Zero
zero = 1.0e-20

## Set the maximum size of xenon
maxSize = 1000001
    
## Create plots
fig = plt.figure()
title = 'Xenon Distribution'
fig.suptitle(title,fontsize=22)
xePlot = plt.subplot(111)

## Create lists of file names to read from, time step number in the file, associated line colors, labels
name = ['/home/sophie/Workspace/xolotl-plsm-build/script/xolotlStop.h5', '/home/sophie/Workspace/xolotl-plsm-build/script/xolotlStop.h5', '/home/sophie/Workspace/xolotl-plsm-build/script/xolotlStop.h5', '/home/sophie/Workspace/xolotl-plsm-build/script/xolotlStop.h5']
timestep = [25,35,45,55]
col = ['black', 'blue', 'magenta', 'green']
lab = ['TS 25', 'TS 35', 'TS 45', 'TS 55']

for i in range(len(name)):
    ## Open the file
    f = h5py.File(name[i], 'r')

    ## Open the concentration group
    groupName ='concentrationsGroup/concentration_0_' + str(timestep[i])
    concGroup = f[groupName]

    ## Read the concentration and index datasets
    concDset = concGroup['concs']
    indexDset = concGroup['concs_startingIndices']

    ## Read the time at the chosen time step
    time = concGroup.attrs['absoluteTime']

    ## Read how many normal and super clusters there are
    networkGroup = f['networkGroup']
    totalSize = networkGroup.attrs['totalSize']

    ## Create the mesh and data array
    x = np.empty([maxSize])
    xeArray = np.empty([maxSize])
    for j in range(maxSize):
        x[j] = j
        xeArray[j] = zero
    
    pos = 0 ## if 0D
    for j in range(indexDset[pos], indexDset[pos+1]):
        ## Skip the moments for now
        if (int(concDset[j][0]) > totalSize - 1): continue
        ## Get the cluster bounds
        groupName = str(concDset[j][0])
        clusterGroup = networkGroup[groupName]
        bounds = clusterGroup.attrs['bounds']
        ## Loop on Xe size
        for l in range(bounds[0], bounds[1]+1):
            ## Fill the array
            xeArray[l] = xeArray[l] + concDset[j][1]

    ## Plot the data
    x = np.delete(x,(0), axis=0)
    xeArray = np.delete(xeArray,(0), axis=0)
    xePlot.plot(x, xeArray, lw=4, color=col[i], label=lab[i], alpha=0.75)

## Some formatting
xePlot.set_xlabel("Cluster Size",fontsize=22)
xePlot.set_ylabel("Concentration (# / nm3)",fontsize=22)
xePlot.set_xlim([1, 1000000])
xePlot.set_ylim([1.0e-16, 1.0e-1])
xePlot.set_xscale('log')
xePlot.set_yscale('log')
xePlot.tick_params(axis='both', which='major', labelsize=20)
## Plot the legends
l2 = xePlot.legend(loc='best')
setp(l2.get_texts(), fontsize=25)

## Show the plots
plt.show()
