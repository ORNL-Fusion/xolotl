#!/usr/bin/env python

import numpy as np
import math
import matplotlib.pyplot as plt
from   pylab import *
import h5py
from matplotlib.colors import LogNorm
#from superMapFile import superMap

## Set the Zero
zero = 1.0e-20

## Set the maximum size of helium/hydrogen
maxSize = 201

## Open the file
f = h5py.File('/home/sophie/Workspace/xolotl-plsm-build/script/xolotlStop.h5', 'r')

## Get the last time step saved in the file
concGroup0 = f['concentrationsGroup']
timestep = concGroup0.attrs['lastTimeStep']
lastLoop = concGroup0.attrs['lastLoop']

## Open the concentration group
groupName ='concentrationsGroup/concentration_' + str(lastLoop) + '_' + str(timestep)
concGroup = f[groupName]

## Read the concentration and index datasets
concDset = concGroup['concs']
indexDset = concGroup['concs_startingIndices']

## Read the time at the chosen time step
time = concGroup.attrs['absoluteTime']

## Read the grid to know which grid point is which depth
gridDset = concGroup['grid']
gridSize = len(gridDset)

## Read how many normal and super clusters there are
networkGroup = f['networkGroup']
totalSize = networkGroup.attrs['totalSize']

## Create the mesh and data array
x = np.empty([maxSize+1, gridSize])
y = np.empty([maxSize+1, gridSize])
heArray = np.empty([maxSize+1, gridSize])
tArray = np.empty([maxSize+1, gridSize])
vArray = np.empty([maxSize+1, gridSize])
for i in range(0, maxSize+1):
    for j in range(0, gridSize):
        x[i][j] = gridDset[j] - gridDset[1]
        y[i][j] = i
        heArray[i][j] = zero
        tArray[i][j] = zero
        vArray[i][j] = zero
        
## Save the composition information for all clusters
superMap = []
for i in range(0, totalSize):
    ## Get the cluster bounds
    groupName = str(i)
    clusterGroup = networkGroup[groupName]
    bounds = clusterGroup.attrs['bounds']
    if (bounds[8] > 0): continue # I case
    temp = [[],[],[]]
    ## Loop on them
    for he in range(bounds[0], bounds[1]+1):
        for d in range(bounds[2], bounds[3]+1):
            for t in range(bounds[4], bounds[5]+1):
                for v in range(bounds[6], bounds[7]+1):
                    ## Look for he size
                    foundSize = False
                    for j in range(0, len(temp[0])):
                        if (he == temp[0][j][0]):
                            ## Add to the weight
                            temp[0][j][1] = temp[0][j][1] + 1
                            foundSize = True
                            break
                    ## Create a new field
                    if (not foundSize):
                        temp[0].append([he, 1])
                    ## Look for t size
                    foundSize = False
                    for j in range(0, len(temp[1])):
                        if (t == temp[1][j][0]):
                            ## Add to the weight
                            temp[1][j][1] = temp[1][j][1] + 1
                            foundSize = True
                            break
                    ## Create a new field
                    if (not foundSize):
                        temp[1].append([t, 1])
                    ## Look for v size
                    foundSize = False
                    for j in range(0, len(temp[2])):
                        if (v == temp[2][j][0]):
                            ## Add to the weight
                            temp[2][j][1] = temp[2][j][1] + 1
                            foundSize = True
                            break
                    ## Create a new field
                    if (not foundSize):
                        temp[2].append([v, 1])
    
    superMap.append(temp)
    
#print('superMap =', superMap)

## Loop on the grid
for j in range(len(indexDset)-1):
    ## Loop on the concentrations
    for i in range(indexDset[j], indexDset[j+1]):
        ## Skip the moments for now
        if (int(concDset[i][0]) > totalSize - 1): continue
        ## Get the concentration
        conc = concDset[i][1]
        ## Loop on each component
        lists = superMap[concDset[i][0]]
        ## He
        for k in range(0, len(lists[0])):
            heSize = lists[0][k][0]
            heArray[heSize][j] = heArray[heSize][j] + (conc * lists[0][k][1])
        ## T
        for k in range(0, len(lists[1])):
            tSize = lists[1][k][0]
            tArray[tSize][j] = tArray[tSize][j] + (conc * lists[1][k][1])
        ## V
        for k in range(0, len(lists[2])):
            vSize = lists[2][k][0]
            vArray[vSize][j] = vArray[vSize][j] + (conc * lists[2][k][1])


## Create plots
fig = plt.figure()
title = 'Helium Concentrations at t = ' + str(time) + ' s'
fig.suptitle(title,fontsize=22)
hePlot = plt.subplot(111)
zMaxPlus = max(heArray.max(), max(tArray.max(), vArray.max()))

## Plot the data
cb1 = hePlot.pcolor(x, y, heArray, norm=LogNorm(vmin=1.0e-10, vmax=zMaxPlus), cmap="Oranges", alpha=1.0)
#cb1 = hePlot.pcolor(x, y, heArray, vmin=-zMaxPlus, vmax=zMaxPlus, cmap="bwr", alpha=1.0)
hePlot.set_xlabel("Depth (nm)",fontsize=22)
hePlot.set_ylabel("Helium Cluster Size",fontsize=22)
hePlot.set_xlim([0.1, gridDset[len(gridDset)-1] - gridDset[1]])
hePlot.set_ylim([0, maxSize])
hePlot.set_xscale('log')
hePlot.tick_params(axis='both', which='major', labelsize=20)

## Make an axis for the colorbar on the right side
cax = fig.add_axes([0.9, 0.1, 0.03, 0.8])
colorbar1 = plt.colorbar(cb1, cax=cax)
colorbar1.set_label("# / nm3",fontsize=22)
colorbar1.ax.tick_params(axis='y',labelsize=20)

#fig.savefig('He_cong.png')

## Create plots
fig2 = plt.figure()
title = 'Tritium Concentrations at t = ' + str(time) + ' s'
fig2.suptitle(title,fontsize=22)
tPlot = plt.subplot(111)

## Plot the data
cb2 = tPlot.pcolor(x, y, tArray, norm=LogNorm(vmin=1.0e-10, vmax=zMaxPlus), cmap="Oranges", alpha=1.0)
#cb2 = tPlot.pcolor(x, y, tArray, vmin=-zMaxPlus, vmax=zMaxPlus, cmap="bwr", alpha=1.0)
tPlot.set_xlabel("Depth (nm)",fontsize=22)
tPlot.set_ylabel("Tritium Cluster Size",fontsize=22)
tPlot.set_xlim([0.1, gridDset[len(gridDset)-1] - gridDset[1]])
tPlot.set_ylim([0, maxSize])
tPlot.set_xscale('log')
tPlot.tick_params(axis='both', which='major', labelsize=20)

## Make an axis for the colorbar on the right side
cax = fig2.add_axes([0.9, 0.1, 0.03, 0.8])
colorbar2 = plt.colorbar(cb2, cax=cax)
colorbar2.set_label("# / nm3",fontsize=22)
colorbar2.ax.tick_params(axis='y',labelsize=20)

#fig2.savefig('T_cong.png')

## Create plots
fig3 = plt.figure()
title = 'Vacancy Concentrations at t = ' + str(time) + ' s'
fig3.suptitle(title,fontsize=22)
vPlot = plt.subplot(111)

## Plot the data
cb3 = vPlot.pcolor(x, y, vArray, norm=LogNorm(vmin=1.0e-10, vmax=zMaxPlus), cmap="Oranges", alpha=1.0)
#cb3 = vPlot.pcolor(x, y, vArray, vmin=-zMaxPlus, vmax=zMaxPlus, cmap="bwr", alpha=1.0)
vPlot.set_xlabel("Depth (nm)",fontsize=22)
vPlot.set_ylabel("Tritium Cluster Size",fontsize=22)
vPlot.set_xlim([0.1, gridDset[len(gridDset)-1] - gridDset[1]])
vPlot.set_ylim([0, maxSize])
vPlot.set_xscale('log')
vPlot.tick_params(axis='both', which='major', labelsize=20)

## Make an axis for the colorbar on the right side
cax = fig3.add_axes([0.9, 0.1, 0.03, 0.8])
colorbar3 = plt.colorbar(cb3, cax=cax)
colorbar3.set_label("# / nm3",fontsize=22)
colorbar3.ax.tick_params(axis='y',labelsize=20)

#fig3.savefig('V_cong.png')

## Show the plots
plt.show()
