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
temperature = 1600.0
k_B = 1.38066e-5

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
maxSize = 501
volume = 4.09e-2 # nm3
maxRadius = int(pow(3.0 * maxSize * volume / (4.0 * math.pi), 1.0/3.0)) + 1


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
    
    minR = pow(3.0 * bounds[2] * volume / (4.0 * math.pi), 1.0/3.0)
    maxR = pow(3.0 * (bounds[3] + 1) * volume / (4.0 * math.pi), 1.0/3.0)
    
    T_r = temperature / 290.0
    avXe = (bounds[0] + bounds[1]) / 2.0
    avVolume = volume * (bounds[2] + bounds[3]) / 2.0
    avR = pow(3.0 * avVolume / (4.0 * math.pi), 1.0/3.0)
    rho_B = avXe / avVolume
    rho_Br = rho_B / 5.06
    delta_R = (-0.00233632 + 0.00117093 * T_r) * rho_Br + (0.00056123 - 0.00015370 * T_r) * rho_Br * rho_Br
    term = avR / (avR + delta_R)
    rho = rho_B * term * term * term
    rho_r = rho / 5.06
    B = 0.261323 - 1.132763 / T_r + 0.028564 / (T_r * T_r)
    C = 0.318315 - 0.038613 / T_r - 0.248067 / (T_r * T_r)
    E = -0.003885 + 0.069736 / T_r + 0.018356 / (T_r * T_r)
    pressure = rho * k_B * temperature * (1.0 + B * rho_r + C * rho_r * rho_r + E * rho_r * rho_r * rho_r * rho_r)
    
    ## Save info for rectangles
    if (concDset[j][1] > 0):
        vertsPos.append(list(zip([minR,maxR,maxR,minR,minR], [bounds[0],bounds[0],bounds[1] + 1,bounds[1] + 1,bounds[0]])))
        cPos = np.append(cPos, max(zero,pressure))
        
        
    
## Create the figure
figPos,axPos=plt.subplots()
axPos.set_ylim(0,maxSize)
axPos.set_xlim(0,maxRadius)

#, norm=matplotlib.colors.LogNorm(vmin=1.0e-16)
## Create the rectangle collection
colPos = PolyCollection(vertsPos, cmap="winter")
colPos.set_array(cPos)
colPos.set_edgecolor('k')
colPos.set_linewidth(0.5)
axPos.add_collection(colPos)

## Set axis labels
axPos.set_xlabel("Bubble Radius [nm]",fontsize=22)
axPos.set_ylabel("Xenon Cluster Size",fontsize=22)
axPos.tick_params(axis='both', which='major', labelsize=20)

## The colorbar
cbarPos = figPos.colorbar(colPos)
cbarPos.set_label("[GPa]",fontsize=22)
cbarPos.ax.tick_params(axis='y',labelsize=20)

plt.show()
