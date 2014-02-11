#!/usr/bin/env python

import numpy as np
import math
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from   pylab import *
from   scipy import stats
from scipy import interpolate
from scipy import optimize
from scipy.stats import norm
 
## Create plot
fig = plt.figure()
HeSpread = plt.subplot(111)

## List to store all the Helium formation energies
HeEnergy = []

## Load helium number, vacancy number, and formation energy from 'formationEnergies.dat'
He, V, formationEnergyComputed = loadtxt('formationEnergies.dat', usecols = (1,0,2) , unpack=True)
 
## Load helium number, vacancy number, and binding energy from 'tungsten.txt'
HeTrue, VTrue, bindingEnergyTrue = loadtxt('../../../benchmarks/tungsten.txt', usecols = (0,1,3) , unpack=True) 

## Loop on possible vacancy numbers
for i in range(1, 50):
    
    ## Filter on the vacancy number
    HeF = He[V == i]
    formationEnergyComputedF = formationEnergyComputed[V == i]
    VF = i 
    
    ## If data exists 
    if len(HeF) > 0:
        
        ## Previous energies
        previousEnergyComputed = formationEnergyComputedF[0]
        
        ## Loop on Helium number
        for k in range(1, len(HeF)):
            
            ## Filter on Helium number and vacancy number 
            bindingEnergyTrueF = bindingEnergyTrue[HeTrue == HeF[k]]
            VTrueF = VTrue[HeTrue == HeF[k]]
            bindingEnergyTrueFF = bindingEnergyTrueF[VTrueF == VF]
            
            ## Compute the Helium formation energy
            energy = float(bindingEnergyTrueFF + formationEnergyComputedF[k] - previousEnergyComputed)
            HeEnergy.append(energy)
            previousEnergyComputed = formationEnergyComputedF[k]
            
           
## Plot the histogram 
data = HeSpread.hist(HeEnergy, 50, color='c', normed=True) 

## Fit a normal distribution to the data
mean, deviation = norm.fit(HeEnergy)

## And plot it
x = np.linspace(min(HeEnergy), max(HeEnergy), 200)
HeSpread.plot(x, norm.pdf(x, mean, deviation), linewidth=5, color='m') 

## Name the plot
title = "Helium Formation Energy: mean = %.2f,  standard deviation = %.2f" % (mean, deviation)
HeSpread.set_title(title, fontsize=20)
        
## Show the plots
plt.show() 
        