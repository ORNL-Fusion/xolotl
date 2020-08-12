#!/usr/bin/env python

import numpy as np
import math
import matplotlib.pyplot as plt
from   pylab import *
import h5py

## Create plots
fig = plt.figure()
plotDist = plt.subplot(111)

## Open the file
f = h5py.File('/home/sophie/Data/Xolotl/PISCES/TRIDYN_0.h5', 'r')
concDset = f['concs']
depth = []
He = []
D = []
T = []
V = []
I = []
## Read the data and put it in lists
for j in range(0, len(concDset)):
    depth.append(concDset[j][0])
    He.append(concDset[j][1]) ## 1 for He, 2 for D, 3 for T, 4 for V, 5 for I
    D.append(concDset[j][2])
    T.append(concDset[j][3])
    V.append(concDset[j][4])
    I.append(concDset[j][5])
plotDist.plot(depth, He, ls='-', lw=4, marker='.', markersize=10, color='orange', alpha=0.5, label="He")
plotDist.plot(depth, D, ls='-', lw=4, marker='.', markersize=10, color='cyan', alpha=0.5, label="D")
plotDist.plot(depth, T, ls='-', lw=4, marker='.', markersize=10, color='green', alpha=0.5, label="T")
plotDist.plot(depth, V, ls='-', lw=4, marker='.', markersize=10, color='blue', alpha=0.5, label="V")
plotDist.plot(depth, I, ls='-', lw=4, marker='.', markersize=10, color='magenta', alpha=0.5, label="I")

## Formatting
plotDist.set_xlabel("Depth [nm]",fontsize=25)
plotDist.set_ylabel("Concentration [atoms/nm3]",fontsize=25)
#plotDist.set_yscale('log')
#plotDist.set_xlim([0.0, 1000.0])
#plotDist.set_ylim([0.0, 0.15])
#plotDist.grid()
plotDist.tick_params(axis='both', which='major', labelsize=25)
plotDist.tick_params(axis='both', which='minor', labelsize=25)

## Plot the legends
l2 = plotDist.legend(loc='best')
setp(l2.get_texts(), fontsize=25)

## Show the plots
plt.show()
