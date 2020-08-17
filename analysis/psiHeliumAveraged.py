#!/usr/bin/env python
import numpy as np
import math
import matplotlib.pyplot as plt
from   pylab import *

## Create plots
fig1 = plt.figure()
conPlot = plt.subplot(111)

## Get the helium data from a retentionOut.txt file
fluence1, He1 = loadtxt('/home/sophie/Data/Xolotl/dynamicGrid/retention_dynamic_no.txt', usecols = (0,1) , unpack=True)

## Give the flux used in the simulation to plot as a function of time
flux = 54000

## Loop on the data
time_limit = 1.0 # period of time in s over which the content is averaged
count = 0
conc = 0.0
previous_time = 0.0
x = []
y = []
for i in range(len(fluence1)):
    time = fluence1[i] / flux
    if (time > time_limit):
        if (count > 0):
            x.append(time_limit)
            y.append(conc * 1.0e-1) # factor to convert from nm-2 to e15 cm-2
            conc = 0.0
            count = 0
        time_limit = time_limit + 1.0
    conc = conc + (He1[i] * (time - previous_time))
    count = count + 1
    previous_time = time
conPlot.plot(x, y, linewidth=2, color='k', alpha=1.0, label='No moving surface')

## Plot the legend
l = conPlot.legend(loc='best', ncol=2)
setp(l.get_texts(), fontsize=25)

## Some shaping
conPlot.set_xlabel("Time (s)",fontsize=25)
conPlot.set_ylabel("Helium content (x $10^{15}$ cm$^{-2}$)",fontsize=25)
#conPlot.set_xlim([0.0, 1000.0])
#conPlot.set_ylim([0.0001,1000])
#conPlot.set_xscale('log')
#conPlot.set_yscale('log')
conPlot.grid()
conPlot.tick_params(axis='both', which='major', labelsize=25)
conPlot.tick_params(axis='both', which='minor', labelsize=25)

## Show the plots
plt.show()
