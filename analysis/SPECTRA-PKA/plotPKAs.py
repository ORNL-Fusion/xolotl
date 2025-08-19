#!/usr/bin/env python

import numpy as np
import math
import matplotlib.pyplot as plt
from   pylab import *
import matplotlib as mpl

## Create plots
fig1 = plt.figure()
pkaPlot = plt.subplot(111)

## Constants
N_A = 6.02214129e23
density = 19.254 #g cm-3
atomic_weight = 183.84

lo, hi, pka = loadtxt('W_PWR_elemental.dat', usecols = (0,1,2) , unpack=True)

x = []
for i in range(len(lo)):
    x.append((lo[i]+hi[i])*5e5)

pkaPlot.plot(x, (pka * N_A * density)/atomic_weight, linewidth=8, ls='-', color='k', alpha=1, label='PWR')

lo, hi, pka = loadtxt('W_fast_higher_elemental.dat', usecols = (0,1,2) , unpack=True)

x = []
for i in range(len(lo)):
    x.append((lo[i]+hi[i])*5e5)

pkaPlot.plot(x, (pka * N_A * density)/atomic_weight, linewidth=6, ls='-', color='magenta', alpha=1, label='fast')

## Plot the legend
l = pkaPlot.legend(loc='best')
setp(l.get_texts(), fontsize=20)

## Some shaping
pkaPlot.set_xlabel("PKAs energy [eV]",fontsize=25)
pkaPlot.set_ylabel("PKAs s$^{-1}$ cm$^{-3}$",fontsize=25)
pkaPlot.set_xlim([1.0, 1.0e6])
pkaPlot.set_ylim([1.0e6, 1.0e13])
pkaPlot.set_xscale('log')
pkaPlot.set_yscale('log')
pkaPlot.tick_params(axis='both', which='major', labelsize=25)
pkaPlot.tick_params(axis='both', which='minor', labelsize=25)

## Show the plots
plt.show()
