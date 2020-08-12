#!/usr/bin/env python

import numpy as np
import math
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.collections import PolyCollection
from   pylab import *
import matplotlib.transforms as mtrans

## Open the Xolotl output file (by default "bubble_#.dat")
amin, amax, bmin, bmax, conc = loadtxt('/home/sophie/Workspace/xolotl-master-build/bubble_211.dat', usecols = (0,1,2,3,4) , unpack=True)

## Init objects
x = np.array([])
y = np.array([])
c = np.array([])
verts = []
zero = 1.0e-20
## Loop on the concentrations to create the rectangle objects
for k in range(len(conc)):
    x = np.append(x, bmin[k])
    y = np.append(y, amin[k])
    awidth = amax[k] - amin[k] + 1
    bwidth = bmax[k] - bmin[k] + 1
    verts.append(list(zip([0,bwidth,bwidth,0,0], [0,0,awidth,awidth,0])))
    c = np.append(c, log(max(zero,conc[k]))) # log scale for the concentrations   
    
## Create the figure
fig,ax=plt.subplots()
plt.ylim(0,amax.max())
plt.xlim(0,bmax.max())

## Process the offsets and create the rectangle collection
offsets = np.c_[x,y]
col = PolyCollection(verts, offsets=offsets, transOffset=mtrans.IdentityTransform(), offset_position="data", cmap="spring")
col.set_array(c)
col.set_edgecolor('k')
col.set_linewidth(0.5)
ax.add_collection(col)

## Set axis labels
ax.set_xlabel("Vacancy Cluster Size",fontsize=22)
ax.set_ylabel("Helium Cluster Size",fontsize=22)
ax.tick_params(axis='both', which='major', labelsize=20)

## The colorbar
cbar = fig.colorbar(col)
cbar.set_label("log(# / nm3)",fontsize=22)
cbar.ax.tick_params(axis='y',labelsize=20)

plt.show()
