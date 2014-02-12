#!/usr/bin/env python
#=======================================================================================
# computeBindingEnergies.py
# Plots three things:
#    - real helium binding energies from the file '../../../benchmarks/tungsten.txt'
#    - computed helium binding energies with real formation energies from the file  
# 'formationEnergies.dat'
#    - computed helium binding energies with fitted formation energies from the file
# 'outputFile.dat' created by formationEnergiesFit.py
# The formula used to compute binding energies is:
#   Eb(He_x, V_y) = Ef(He_x-1, V_y) + Ef(He) - Ef(He_x, V_y)
# where you have to set Ef(He)
# Can compute the squared distance between real data and computed data by uncommenting 
# the end of this script
#=======================================================================================

import numpy as np
import math
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from   pylab import *
from   scipy import stats
from scipy import interpolate

## Set the value of the Helium formation energy to use
HeFormation = 6.15
 
## Create plots
fig = plt.figure()
energiesComputed = plt.subplot(312)
energiesTrue = plt.subplot(311)
energiesFit = plt.subplot(313)

## Load helium number, vacancy number, and formation energy from 'formationEnergies.dat'
He, V, formationEnergy = loadtxt('formationEnergies.dat', usecols = (1,0,2) , unpack=True)

## List to store which vacancy numbers are used
VList = []

## Loop on possible vacancy numbers
for i in range(1, 50):
    
    ## Filter on the vacancy number
    HeFiltered = He[V == i]
    formationEnergyFiltered = formationEnergy[V == i]
    VFiltered = i 
   
    ## If data exists 
    if len(HeFiltered) > 0:

        ## Initialize the previous value of formation energy
        energyPrevious = formationEnergyFiltered[0]

        ## Declare lists where binding energies will be stored
        bindingEnergy = []  
        
        ## Loop on all the HeFiltered elements
        for k in range(1, len(HeFiltered)):
            
            ## Compute binding energy assuming Ef(He, 0) = HeFormation
            value = energyPrevious + HeFormation - formationEnergyFiltered[k]
            bindingEnergy.append(value)
            energyPrevious = formationEnergyFiltered[k]
            
        ## Remove first item of the list because the binding energy cannot be computed for it
        HeFiltered = np.delete(HeFiltered, 0, axis=0)
        VFiltered = np.delete(VFiltered, 0, axis=0)
        
        ## Plot binding energy = f(He/V)
        energiesComputed.plot(HeFiltered/VFiltered, bindingEnergy, color=(0.02*i,0,1-0.02*i), linewidth=1.5)
        
        ## Store which vacancy numbers were used
        VList.append(VFiltered)

## Load helium number, vacancy number, and formation energy from 'outputFile.dat'
He, V, formationEnergy = loadtxt('outputFile.dat', usecols = (1,0,2) , unpack=True)

## Loop on possible vacancy numbers
for i in range(1, 50):
    
    ## Filter on the vacancy number
    HeFiltered = He[V == i]
    formationEnergyFiltered = formationEnergy[V == i]
    VFiltered = i 
   
    ## If data exists 
    if len(HeFiltered) > 0:

        ## Initialize the previous value of formation energy
        energyPrevious = formationEnergyFiltered[0]

        ## Declare lists where binding energies will be stored
        bindingEnergy = []  
        
        ## Loop on all the HeFiltered elements
        for k in range(1, len(HeFiltered)):
            
            ## Compute binding energy assuming Ef(He, 0) = HeFormation
            value = energyPrevious + HeFormation - formationEnergyFiltered[k]
            bindingEnergy.append(value)
            energyPrevious = formationEnergyFiltered[k]
            
        ## Remove first item of the list because the binding energy cannot be computed for it
        HeFiltered = np.delete(HeFiltered, 0, axis=0)
        VFiltered = np.delete(VFiltered, 0, axis=0)
        
        ## Plot binding energy = f(He/V)
        energiesFit.plot(HeFiltered/VFiltered, bindingEnergy, color=(0.02*i,0,1-0.02*i), linewidth=1.5)
 
## Load helium number, vacancy number, and binding energy from 'tungsten.txt'
He, V, bindingEnergy = loadtxt('../../../benchmarks/tungsten.txt', usecols = (0,1,3) , unpack=True)  

## Loop on possible vacancy numbers
for i in range(0, len(VList)):
    
    ## Filter on the vacancy number
    HeFiltered = He[V == VList[i]]
    bindingEnergyFiltered = bindingEnergy[V == VList[i]]
    VFiltered = VList[i]
   
    ## If data exists 
    if len(HeFiltered) > 0:
        
        ## Plot binding energy = f(He/V)
        energiesTrue.plot(HeFiltered/VFiltered, bindingEnergyFiltered, color=(0.02*VList[i],0,1-0.02*VList[i]), linewidth=1.5)       

## Set the same X and Y scale for all plots
energiesTrue.set_xlim([0, 8])
energiesTrue.set_ylim([1, 7])
energiesComputed.set_xlim([0, 8])
energiesComputed.set_ylim([1, 7])
energiesFit.set_xlim([0, 8])
energiesFit.set_ylim([1, 7])

## Title and axis
energiesTrue.set_title("Energies from Benchmark", fontsize=20)
energiesTrue.set_ylabel("Helium binding energy",fontsize=16)
energiesTrue.grid()
energiesComputed.set_title("Energies computed from formation energies", fontsize=20)
energiesComputed.set_ylabel("Helium binding energy",fontsize=16)
energiesComputed.grid()
energiesFit.set_title("Energies computed from fitted formation energies", fontsize=20)
energiesFit.set_xlabel("Helium/Vacancy number",fontsize=16)
energiesFit.set_ylabel("Helium binding energy",fontsize=16)
energiesFit.grid()

## Uncomment the following to compute the least squares
#   
# ## Compute Least Squares
# leastSquareComputed = 0.
# leastSquareFit = 0.
#   
# ## Load helium number, vacancy number, and formation energy from 'formationEnergies.dat'
# He, V, formationEnergyComputed = loadtxt('formationEnergies.dat', usecols = (1,0,2) , unpack=True)
#   
# ## Load formation energy from 'outputFile.dat'
# He, V, formationEnergyFit = loadtxt('outputFile.dat', usecols = (1,0,2) , unpack=True)
#    
# ## Load helium number, vacancy number, and binding energy from 'tungsten.txt'
# HeTrue, VTrue, bindingEnergyTrue = loadtxt('../../../benchmarks/tungsten.txt', usecols = (0,1,3) , unpack=True) 
#   
# ## Loop on possible vacancy numbers
# for i in range(1, 50):
#       
#     ## Filter on the vacancy number
#     HeF = He[V == i]
#     formationEnergyComputedF = formationEnergyComputed[V == i]
#     formationEnergyFitF = formationEnergyFit[V == i]
#     VF = i 
#       
#     ## If data exists 
#     if len(HeF) > 0:
#           
#         ## Previous energies
#         previousEnergyComputed = formationEnergyComputedF[0]
#         previousEnergyFit = formationEnergyFitF[0]
#           
#         ## Loop on Helium number
#         for k in range(1, len(HeF)):
#               
#             ## Filter on Helium number and vacancy number 
#             bindingEnergyTrueF = bindingEnergyTrue[HeTrue == HeF[k]]
#             VTrueF = VTrue[HeTrue == HeF[k]]
#             bindingEnergyTrueFF = bindingEnergyTrueF[VTrueF == VF]
#               
#             ## Compute binding energies
#             bindingEnergyComputed = previousEnergyComputed + HeFormation - formationEnergyComputedF[k]
#             previousEnergyComputed = formationEnergyComputedF[k]
#             bindingEnergyFit = previousEnergyFit + HeFormation - formationEnergyFitF[k]
#             previousEnergyFit = formationEnergyFitF[k]
#               
#             ## Compute the squared distance
#             rComputed = math.pow((bindingEnergyComputed - bindingEnergyTrueFF), 2)
#             rFit = math.pow((bindingEnergyFit - bindingEnergyTrueFF), 2)
#               
#             ## Add to the total value
#             leastSquareComputed = leastSquareComputed + rComputed
#             leastSquareFit = leastSquareFit + rFit
#   
# ## Print the least squares
# print "Distance from real data to computed ones: ", leastSquareComputed
# print "Distance from real data to fitted ones: ", leastSquareFit
        
## Show the plots
plt.show() 
        