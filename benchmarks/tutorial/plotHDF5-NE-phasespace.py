#!/usr/bin/env python

import numpy as np
import math
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.collections import PolyCollection
from   pylab import *
import h5py
import matplotlib.transforms as mtrans

## Open the file
f = h5py.File('/home/sophie/Workspace/xolotl-NE-build/fullPhasespace/xolotlStop.h5', 'r')

## Get the TS
concGroup = f['concentrationsGroup']
timestep = concGroup.attrs['lastTimeStep']
lastLoop = concGroup.attrs['lastLoop']
#timestep = 415

## Open the concentration group
groupName ='concentrationsGroup/concentration_' + str(lastLoop) + '_' + str(timestep)
concGroup = f[groupName]

## Read the concentration and index datasets
concDset = concGroup['concs']
indexDset = concGroup['concs_startingIndices']

## Read the time at the chosen time step
time = concGroup.attrs['absoluteTime']

print(time)

## Read how many normal and super clusters there are
networkGroup = f['networkGroup']
totalSize = networkGroup.attrs['totalSize']

## Init objects
cPos, cNeg = np.array([]), np.array([])
vertsPos, vertsNeg = [], []
zero = 1.0e-20
maxSize = 51

## Create the data arrays
x = np.empty([maxSize+1])
xeArray = np.empty([maxSize+1])
vArray = np.empty([maxSize+1])
for i in range(0, maxSize+1):
    x[i] = i
    xeArray[i] = zero
    vArray[i] = zero

## Loop on the clusters
pos = 0
for j in range(indexDset[pos], indexDset[pos+1]):
    ## Skip the moments for now
    if (int(concDset[j][0]) > totalSize - 1): continue
    ## Get the cluster bounds
    groupName = str(concDset[j][0])
    clusterGroup = networkGroup[groupName]
    bounds = clusterGroup.attrs['bounds']
    ## Skip I
    if (bounds[4] > 0): continue
    
    ## Save info for rectangles
    if (concDset[j][1] > 0):
        vertsPos.append(list(zip([bounds[2],bounds[3] + 1,bounds[3] + 1,bounds[2],bounds[2]], [bounds[0],bounds[0],bounds[1] + 1,bounds[1] + 1,bounds[0]])))
        cPos = np.append(cPos, max(zero,concDset[j][1]))
        
        ## Loop on V size
        for k in range(bounds[2], bounds[3]+1):
            vArray[k] = vArray[k] + concDset[j][1] * k * (bounds[1]+1 - bounds[0])
        ## Loop on Xe size
        for l in range(bounds[0], bounds[1]+1):
            ## Fill the arrays
            xeArray[l] = xeArray[l] + concDset[j][1] * l * (bounds[3]+1 - bounds[2])
        
    else:
        vertsNeg.append(list(zip([bounds[2],bounds[3] + 1,bounds[3] + 1,bounds[2],bounds[2]], [bounds[0],bounds[0],bounds[1] + 1,bounds[1] + 1,bounds[0]])))
        cNeg = np.append(cNeg, max(zero,-concDset[j][1]))
        
    
## Create the figure
figPos,axPos=plt.subplots()
axPos.set_ylim(0,maxSize)
axPos.set_xlim(0,maxSize)
figNeg,axNeg=plt.subplots()
axNeg.set_ylim(0,maxSize)
axNeg.set_xlim(0,maxSize)
fig = plt.figure()
intPlot = plt.subplot(111)

## Create the rectangle collection
colPos = PolyCollection(vertsPos, cmap="spring", norm=matplotlib.colors.LogNorm(vmin=1.0e-16))
colPos.set_array(cPos)
colPos.set_edgecolor('k')
colPos.set_linewidth(0.5)
axPos.add_collection(colPos)
colNeg = PolyCollection(vertsNeg, cmap="winter", norm=matplotlib.colors.LogNorm(vmin=1.0e-16))
colNeg.set_array(cNeg)
colNeg.set_edgecolor('k')
colNeg.set_linewidth(0.5)
axNeg.add_collection(colNeg)
intPlot.plot(x, xeArray, linewidth=5, color="blue", label="Xe Size")
intPlot.plot(x, vArray, linewidth=5, color="orange", label="V Size")

## Set axis labels
axPos.set_xlabel("Vacancy Cluster Size",fontsize=22)
axPos.set_ylabel("Xenon Cluster Size",fontsize=22)
axPos.tick_params(axis='both', which='major', labelsize=20)
axNeg.set_xlabel("Vacancy Cluster Size",fontsize=22)
axNeg.set_ylabel("Xenon Cluster Size",fontsize=22)
axNeg.tick_params(axis='both', which='major', labelsize=20)

intPlot.set_xlabel("Cluster Size",fontsize=22)
intPlot.set_ylabel("Concentration [# / nm3]",fontsize=22)
intPlot.set_xlim([1, maxSize])
intPlot.set_ylim([1.0e-16, 1.0e-3])
#intPlot.set_xscale('log')
intPlot.set_yscale('log')
intPlot.tick_params(axis='both', which='major', labelsize=20)
## Plot the legend
l = intPlot.legend(loc='best')
setp(l.get_texts(), fontsize=16)

## The colorbar
cbarPos = figPos.colorbar(colPos)
cbarPos.set_label("[# / nm3]",fontsize=22)
cbarPos.ax.tick_params(axis='y',labelsize=20)
if (len(vertsNeg)> 0):
    cbarNeg = figNeg.colorbar(colNeg)
    cbarNeg.set_label("[# / nm3]",fontsize=22)
    cbarNeg.ax.tick_params(axis='y',labelsize=20)

plt.show()
