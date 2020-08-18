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
name = ['/home/sophie/Workspace/xolotl-reso-build/baker/network_2073K.h5', '/home/sophie/Workspace/xolotl-reso-build/baker/network_2073K_full.h5', '/home/sophie/Workspace/xolotl-reso-build/baker/network_2073K_one.h5', '/home/sophie/Workspace/xolotl-reso-build/baker/network_2073K_partial.h5']
timestep = [436,511,436,436]
col = ['black', 'blue', 'magenta', 'green', 'red']
lab = ['2073 K, standard', '2073 K, Turnbull re-solution', '2073 K, Setyawan re-solution (one)', '2073 K, Setyawan re-solution (partial)', 'test']

for i in range(len(name)):
    ## Open the file
    f = h5py.File(name[i], 'r')

    ## Open the concentration group
    groupName ='concentrationsGroup/concentration_' + str(timestep[i])
    concGroup = f[groupName]

    ## Read the concentration and index datasets
    concDset = concGroup['concs']
    indexDset = concGroup['concs_startingIndices']

    ## Read the time at the chosen time step
    time = concGroup.attrs['absoluteTime']

    ## Read the composition index to know which cluster is what
    compDset = f['headerGroup/composition']

    ## Read how many normal and super clusters there are
    networkGroup = f['networkGroup']
    normalSize = networkGroup.attrs['normalSize']
    superSize = networkGroup.attrs['superSize']

    ## Create the mesh and data array
    x = np.empty([maxSize])
    xeArray = np.empty([maxSize])
    for j in range(maxSize):
        x[j] = j
        xeArray[j] = zero
    
    pos = 0 ## if 0D
    for j in range(indexDset[pos], indexDset[pos+1]):
        ## Skip the moments for now
        if (int(concDset[j][0]) > len(compDset) - 1): continue
        ## Normal clusters
        if (int(concDset[j][0]) < normalSize):
            ## Get the helium and hydrogen sizes of this cluster
            xeSize = compDset[int(concDset[j][0])][0]
            ## Fill the arrays
            if (concDset[j][1] > 0.0): xeArray[xeSize] = xeArray[xeSize] + concDset[j][1]
        ## Super clusters
        else:
            ## Loop on the number of clusters it contains
            clusterGroup = networkGroup[str(concDset[j][0])]
            nTot = clusterGroup.attrs['nTot']
            numXe = clusterGroup.attrs['numAtom']
            if (concDset[j][1] > 0.0):
                for l in range(int(numXe - nTot)+1, int(numXe)+1):
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
