#!/usr/bin/env python
#=======================================================================================
# xenonPlotting.py
# Plots the xenon retention values obtained with Xolotl
#=======================================================================================

import numpy as np
import math
import matplotlib.pyplot as plt
from   pylab import *

## Create plots
fig1 = plt.figure()
densityPlot = plt.subplot(111)

time1, I1, V1, V2, Xe1 = loadtxt('/home/sophie/Workspace/xolotl-NE-build/fullPhasespace/retentionOut.txt', usecols = (0,2,3,4,1) , unpack=True)

densityPlot.plot(time1, I1, lw=4, color='orange', ls='-', label='I1', alpha=0.7, marker='o')
densityPlot.plot(time1, V1, lw=4, color='blue', ls='-', label='V1', alpha=0.7, marker='o')
densityPlot.plot(time1, V2, lw=4, color='magenta', ls='-', label='V2', alpha=0.7, marker='o')
densityPlot.plot(time1, Xe1, lw=4, color='black', ls='-', label='Xe', alpha=0.7, marker='o')


## Plot the legend
l = densityPlot.legend(loc="best")
setp(l.get_texts(), fontsize=20)

## Some shaping
densityPlot.set_xlabel("Time (s)",fontsize=25)
densityPlot.set_ylabel("Concentration (nm$^{-3}$)",fontsize=25)
densityPlot.set_xscale('log')
densityPlot.set_yscale('log')
densityPlot.grid()
densityPlot.tick_params(axis='both', which='major', labelsize=25)
densityPlot.tick_params(axis='both', which='minor', labelsize=25)

## Show the plots
plt.show()
