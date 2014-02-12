#!/usr/bin/env python
#=======================================================================================
# confidenceInterval.py
# Computes confidence intervals on the helium binding energy for each vacancy number by:
#    - calculating the distance between real (from the file 
# '../../../benchmarks/tungsten.txt') and fitted data (computed from the file 
# 'outputFile.dat' created by formationEnergiesFit.py)
#    - plotting then in a histogram
#    - fitting this histogram with a normal distribution to obtain the mean and the 
# standard deviation (in a perfect world this mean should be 0)
# Plots both the real binding energies and the intervals obtained to compare
#=======================================================================================

import numpy as np
import math
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from   pylab import *
from   scipy import stats
from scipy import interpolate
from scipy.stats import norm

## Set the value of the Helium formation energy to use
HeFormation = 6.15

## Create plots
fig1 = plt.figure()
gaussianExample = plt.subplot(121)
gaussianDeviation = plt.subplot(222)
gaussianMean = plt.subplot(224)
 
## Create plots
fig2 = plt.figure()
# energiesTrue = plt.subplot(211)
energies = plt.subplot(111)

## List to store all the standard deviations and means
deviationList = []
meanList = []

## List to store which vacancy numbers are used
VList = []
  
## Load formation energy from 'outputFile.dat'
He, V, formationEnergyFit = loadtxt('outputFile.dat', usecols = (1,0,2) , unpack=True)
   
## Load helium number, vacancy number, and binding energy from 'tungsten.txt'
HeTrue, VTrue, bindingEnergyTrue = loadtxt('../../../benchmarks/tungsten.txt', usecols = (0,1,3) , unpack=True) 
  
## Loop on possible vacancy numbers
for i in range(1, 50):
      
    ## Filter on the vacancy number
    HeF = He[V == i]
    formationEnergyFitF = formationEnergyFit[V == i]
    VF = i 
    bindingEnergyTrueF = bindingEnergyTrue[VTrue == i]
    HeTrueF = HeTrue[VTrue == i]
    VTrueF = i
    
    ## Declare the list of distance between the fit and real data
    distance = []
      
    ## If data exists 
    if len(HeF) > 0:
          
        ## Previous energies
        previousEnergyFit = formationEnergyFitF[0]

        ## Declare lists where binding energies will be stored
        bindingFit = []  
          
        ## Loop on Helium number
        for k in range(1, len(HeF)):
              
            ## Filter on Helium number and vacancy number 
            bindingEnergyTrueFF = bindingEnergyTrueF[HeTrueF == HeF[k]]
              
            ## Compute binding energies
            bindingEnergyFit = previousEnergyFit + HeFormation - formationEnergyFitF[k]
            previousEnergyFit = formationEnergyFitF[k]
            bindingFit.append(bindingEnergyFit)
            
            ## Add the distance to the list
            distance.append(float(bindingEnergyTrueFF - bindingEnergyFit))
            
        ## Plot the distance as an histogram
        gaussianExample.hist(distance, 10, color='c', alpha=0.2, normed=True) 
        
        ## Fit a normal distribution to the data
        mean, deviation = norm.fit(distance)

        ## And plot it
        x = np.linspace(min(distance), max(distance), 200)
        gaussianExample.plot(x, norm.pdf(x, mean, deviation), linewidth=3, color='m', alpha=0.3) 
        
        ## Plot binding energy = f(He/V)
        energies.plot(HeTrueF/VTrueF, bindingEnergyTrueF, color=(0.02*i,0,1-0.02*i), linewidth=1.5) 
            
        ## Remove first item of the list because the binding energy cannot be computed for it
        HeF = np.delete(HeF, 0, axis=0)
        
        ## Add confidence interval
        lowerBound = np.sum([bindingFit, [mean-deviation]*len(bindingFit)], axis=0)
        higherBound = np.sum([bindingFit, [mean+deviation]*len(bindingFit)], axis=0)
        
        ## Plot confidence intervals
        energies.fill_between(HeF/VF, lowerBound, higherBound, alpha=0.25, facecolor=(0.02*i,0,1-0.02*i), linewidth=0.)
        
        ## Scatter plots
        gaussianMean.scatter(VF, mean, color="purple", s=100)
        gaussianDeviation.scatter(VF, deviation, color="coral", s=100)
     

## To have understandable plots
gaussianExample.set_title("Distance between real and fitted data", fontsize=20)
gaussianDeviation.set_title("Standard Deviation", fontsize=20)
gaussianDeviation.set_xlabel("Vacancy number",fontsize=16)
gaussianDeviation.set_xlim([0, 50])
gaussianDeviation.grid()
gaussianMean.set_title("Mean", fontsize=20)
gaussianMean.set_xlabel("Vacancy number",fontsize=16)
gaussianMean.set_xlim([0, 50])
gaussianMean.grid()

## Set the same X and Y scale for all plots
energies.set_xlim([0, 8])
energies.set_ylim([1, 7])

## Title and axis
energies.set_title("Energies computed from fitted formation energies", fontsize=20)
energies.set_xlabel("Helium/Vacancy number",fontsize=16)
energies.set_ylabel("Helium binding energy",fontsize=16)
      
## Show the plots
plt.show() 
        