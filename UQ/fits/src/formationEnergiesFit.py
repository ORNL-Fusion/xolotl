#!/usr/bin/env python
#=======================================================================================
# formationEnergiesFit.py
# Takes formation energies from the file  'formationEnergies.dat' and fits them with 
# polynomials in two steps:
#    - fits formation energies = f(He/V) for each vacancy number with a piece with 
# function (you can choose each order of polynomials and the position of the separation)
#    - then fits the obtained parameters = f(V) with polynomials (where you can also 
# choose the order)
# Defines the 2D resultFunc corresponding to the fits and write its results in 
# 'outputFile.dat' which has exactly the same format as 'formationEnergies.dat'
#=======================================================================================

import numpy as np
import math
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from   pylab import *
from   scipy import stats
from scipy import interpolate
import numpy.polynomial.polynomial as poly

## Decide the polynomial order for 'helium' fit and 'vacancy' fit and the separation
separation = .9
HeFitOrderLow = 2
VFitOrderLow = 3
HeFitOrderHigh = 3
VFitOrderHigh = 3

## Load helium number, vacancy number, and formation energy from 'formationEnergies.dat'
He, V, formationEnergy = loadtxt('formationEnergies.dat', usecols = (1,0,2) , unpack=True) 

## Declare a list of color to have plots looking better 
colors = ['r', 'b', 'c', 'k', 'g', 'm', 'y']

## Create plots
fig = plt.figure()
energies = plt.subplot(121)
paramsLow = plt.subplot(222)
paramsHigh = plt.subplot(224)
# results = plt.subplot(212, projection='3d')

## Declare lists where all fit parameters will be stored
paramListLow = []
paramListHigh = []
paramFinalLow = []
paramFinalHigh = []

## List to store which vacancy numbers are used
VList = []

## Loop on possible vacancy numbers
for i in range(1, 50):
    
    ## Filter on the vacancy number
    HeFiltered = He[V == i]
    formationEnergyFiltered = formationEnergy[V == i]
    VFiltered = i 
    
    ## Separating data at He/V = separation
    HeLow= HeFiltered[HeFiltered/VFiltered <= separation]
    formationEnergyLow = formationEnergyFiltered[HeFiltered/VFiltered <= separation]
    HeHigh= HeFiltered[HeFiltered/VFiltered > separation]
    formationEnergyHigh = formationEnergyFiltered[HeFiltered/VFiltered > separation]
   
    ## If data exists 
    if len(HeFiltered) > 0:
        
        ## Fit filtered formationEnergy = f(He/V) with a polynomial of order HeFitOrderLow for the low part
        fitLow = poly.polyfit(HeLow/VFiltered, formationEnergyLow, HeFitOrderLow)
        ## Get the fit function
        fitFuncLow = poly.Polynomial(fitLow)
        
        ## Fit filtered formationEnergy = f(He/V) with a polynomial of order HeFitOrderHigh for the high part
        fitHigh = poly.polyfit(HeHigh/VFiltered, formationEnergyHigh, HeFitOrderHigh)
        ## Get the fit function
        fitFuncHigh = poly.Polynomial(fitHigh)
        
        ## Plot both data and the fit
        energies.plot(HeLow/VFiltered, fitFuncLow(HeLow/VFiltered), 'r--', color=(0.02*i,0,1), linewidth=3.)
        energies.plot(HeHigh/VFiltered, fitFuncHigh(HeHigh/VFiltered), 'r--', color=(0.02*i,0,1), linewidth=3.)
        energies.plot(HeFiltered/VFiltered, formationEnergyFiltered, color=(0.02*i,0,1), linewidth=2.)
        
        ## Store fit parameters
        paramListLow.append(fitLow)
        paramListHigh.append(fitHigh)
        VList.append(VFiltered)
        
        ## Loop on the order of the polynomial
        for k in range(0, len(fitLow)):
            
            ## Plot the fit parameter of order k = f(V)
            paramsLow.scatter(VFiltered, fitLow[k], color=colors[6-k], s=100, alpha=0.7)
        
        ## Loop on the order of the polynomial
        for k in range(0, len(fitHigh)):
            
            ## Plot the fit parameter of order k = f(V)
            paramsHigh.scatter(VFiltered, fitHigh[k], color=colors[k], s=100, alpha=0.7)

## Loop on the order of the polynomial for the Low part
for i in range(0, len(paramListLow[0])):
    
    ## Fit the parameters of the same order for different vacancy numbers
    ## with a polynomial of order VFitOrderLow
    fit= poly.polyfit(VList, [row[i] for row in paramListLow], VFitOrderLow)
    ## Get the new fit function
    fitFunc = poly.Polynomial(fit)
    
    ## Plot the fit and store the parameters
    paramsLow.plot(range(1,50), fitFunc(range(1,50)), color=colors[6-i])
    paramFinalLow.append(fit)

## Loop on the order of the polynomial for the High part
for i in range(0, len(paramListHigh[0])):
    
    ## Fit the parameters of the same order for different vacancy numbers
    ## with a polynomial of order VFitOrderHigh
    fit = poly.polyfit(VList, [row[i] for row in paramListHigh], VFitOrderHigh)
    ## Get the new fit function
    fitFunc = poly.Polynomial(fit)
    
    ## Plot the fit and store the parameters
    paramsHigh.plot(range(1,50), fitFunc(range(1,50)), color=colors[i])
    paramFinalHigh.append(fit)

## Definition of the obtained 2D fitted function
def resultFunc(x, y):
    f = 0
    
    ## Piecewise function
    if x <= separation :
        for i in range(0, len(paramListLow[0])):
            for k in range(0, len(paramFinalLow[0])):
                f += paramFinalLow[i][k] * math.pow(x,i) * math.pow(y,k)
            
    else :          
        for i in range(0, len(paramListHigh[0])):
            for k in range(0, len(paramFinalHigh[0])):
                f += paramFinalHigh[i][k] * math.pow(x,i) * math.pow(y,k)
            
    return f

## Open 'outputFile.dat' where results of the fitted function will be printed
outputFile = open('outputFile.dat', 'w')

## Loop on all the elements in 'HeVexpEHe.dat'
for i in range(0, len(He)):
    
    ## Compute the new formation energies
    value = resultFunc(He[i]/V[i], V[i])
        
    ## Write in the output file
    outputFile.write("%d %d %s \n" %(V[i], He[i], value))
    
## Close the output file
outputFile.close()

## Load helium number, vacancy number, and formation energy from the newly created file
He, V, formationEnergy = loadtxt('outputFile.dat', usecols = (1,0,2) , unpack=True) 

## Loop on possible vacancy numbers
for i in range(1, 50):
    
    ## Filter on the vacancy number
    HeFiltered = He[V == i]
    formationEnergyFiltered = formationEnergy[V == i]
    VFiltered = i 
   
    ## If data exists 
    if len(HeFiltered) > 0:
        
        ## Plot both data and the fit
        energies.plot(HeFiltered/VFiltered, formationEnergyFiltered, color=(0.02*i,0,1), linewidth=1.)

## To have understandable plots
paramsLow.set_title("Fit parameters for He/V <= %.1f" % (separation), fontsize=20)
paramsLow.set_xlabel("Vacancy number",fontsize=16)
paramsLow.set_xlim([0, 50])
paramsHigh.set_title("Fit parameters for He/V > %.1f" % (separation), fontsize=20)
paramsHigh.set_xlabel("Vacancy number",fontsize=16)
paramsHigh.set_xlim([0, 50])
energies.set_title("Formation Energies", fontsize=20)
energies.set_xlabel("Helium/Vacancy number",fontsize=16)
energies.set_xlim([0, 6])
energies.set_ylim([0, 100])

## Show the plots
plt.show()