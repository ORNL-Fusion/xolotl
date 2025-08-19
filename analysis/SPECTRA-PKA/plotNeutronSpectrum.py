#!/usr/bin/env python

import numpy as np
import math
import matplotlib.pyplot as plt
from   pylab import *
import matplotlib as mpl

## Create plots
fig1 = plt.figure()
pkaPlot = plt.subplot(111)

all_data = loadtxt('fluxes_specter.dat', skiprows = 3, unpack=True)
N = len(all_data) // 2
x = all_data[:N]
y = all_data[N+1:]

#for i in range(len(y)):
#    if (i < 450): print("0.0")
#    else: print(y[i])

pkaPlot.plot(x, y, linewidth=8, ls='-', color='k', alpha=1)

## Some shaping
pkaPlot.set_xlabel("Neutron energy [MeV]",fontsize=25)
pkaPlot.set_ylabel("Flux [A.U.]",fontsize=25)
#pkaPlot.set_xlim([1.0, 1.0e6])
#pkaPlot.set_ylim([1.0e6, 1.0e13])
pkaPlot.set_xscale('log')
#pkaPlot.set_yscale('log')
pkaPlot.tick_params(axis='both', which='major', labelsize=25)
pkaPlot.tick_params(axis='both', which='minor', labelsize=25)

## Show the plots
plt.show()
