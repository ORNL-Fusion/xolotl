#!/usr/bin/env python
#=======================================================================================
# heliumRetention.py
# Plots the helium retention values as a function of the fluence obtained with Xolotl
#=======================================================================================

import numpy as np
import math
import matplotlib.pyplot as plt
from   pylab import *

## Create plots
fig1 = plt.figure()
conPlot = plt.subplot(111)

## Load the data
fluence1, retention1 = loadtxt('/path/to/data/retention_1.dat', usecols = (0,1) , unpack=True)
fluence2, retention2 = loadtxt('/path/to/data/retention_2.dat', usecols = (0,1) , unpack=True)

## Fill the plot with data
conPlot.scatter(fluence1, retention1, s=100, color='k', label='1')
conPlot.scatter(fluence2, retention2, s=100, color='b', label='1')

## Plot the legend
l = conPlot.legend(loc='best')
setp(l.get_texts(), fontsize=20)
 
## Some cometics
# conPlot.set_title("Title",fontsize=30)
conPlot.set_xlabel("Fluence (x $10^{18}$ He m$^{-2}$)",fontsize=25)
conPlot.set_ylabel("Helium retention (%)",fontsize=25)
conPlot.set_xlim([0.0, 100.0])
conPlot.set_ylim([0.0, 40.0])
# conPlot.set_yscale('log') ## For a log scale
conPlot.grid()
conPlot.tick_params(axis='both', which='major', labelsize=25)
conPlot.tick_params(axis='both', which='minor', labelsize=25)
 
## Show the plots
plt.show()