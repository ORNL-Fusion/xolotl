#!/usr/bin/env python
#=======================================================================================
# heliumCumul.py
# Plots the helium cumulative distribution values as a function of the depth obtained with Xolotl
#=======================================================================================

import numpy as np
import math
import matplotlib.pyplot as plt
from   pylab import *

## Create plot
fig1 = plt.figure()
cumulPlot = plt.subplot(111)

## Load the data
depth1, cumul1 = loadtxt('/path/to/data/heliumCumul_1.dat', usecols = (0,1) , unpack=True)
depth2, cumul2 = loadtxt('/path/to/data/heliumCumul_2.dat', usecols = (0,1) , unpack=True)

## Fill the plot with data
cumulPlot.scatter(depth1, cumul1 / (cumul1[len(cumul1) - 1]), s=100, color='k', label='1')
cumulPlot.scatter(depth2, cumul2 / (cumul2[len(cumul2) - 1]), s=100, color='b', label='2')

## Plot the legend
l = conPlot.legend(loc=4)
setp(l.get_texts(), fontsize=20)
 
## Some cosmetics
#cumulPlot.set_title("Title",fontsize=30)
cumulPlot.set_xlabel("depth (nm)",fontsize=25)
cumulPlot.set_ylabel("Cumulative Fraction of Helium",fontsize=25)
cumulPlot.set_xlim([0.0, 25.0])
cumulPlot.set_ylim([0.0, 1.0])
cumulPlot.grid()
cumulPlot.tick_params(axis='both', which='major', labelsize=25)
cumulPlot.tick_params(axis='both', which='minor', labelsize=25)
 
## Show the plot
plt.show()