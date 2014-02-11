#!/usr/bin/env python

import numpy as np
import math
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from   pylab import *
from   scipy import stats
from scipy import interpolate
import numpy.polynomial.polynomial as poly

## Decide the polynomial order for 'helium' fit and 'vacancy' fit
HeFitOrder = 6
VFitOrder = 4

## Load helium number, vacancy number, and formation energy from 'formationEnergies.dat'
He, V, formationEnergy = loadtxt('formationEnergies.dat', usecols = (1,0,2) , unpack=True) 

## Declare a list of color to have plots looking better 
colors = ['r', 'g', 'c', 'k', 'b', 'm', 'y']

## Create plots
fig = plt.figure()
energies = plt.subplot(121)
params = plt.subplot(122)
# results = plt.subplot(212, projection='3d')

## Declare lists where all fit parameters will be stored
paramList = []
paramFinalList = []

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
        
        ## Fit filtered formationEnergy = f(He/V) with a polynomial of order HeFitOrder
        fit = poly.polyfit(HeFiltered/VFiltered, formationEnergyFiltered, HeFitOrder)
        ## Get the fit function
        fitFunc = poly.Polynomial(fit)
        
        ## Plot both data and the fit
        energies.plot(HeFiltered/VFiltered, fitFunc(HeFiltered/VFiltered), 'r--', color=(0.02*i,0,1), linewidth=3.)
        energies.plot(HeFiltered/VFiltered, formationEnergyFiltered, color=(0.02*i,0,1), linewidth=2.)
        
        ## Store fit parameters
        paramList.append(fit)
        VList.append(VFiltered)
        
        ## Loop on the order of the polynomial
        for k in range(0, len(fit)):
            
            ## Plot the fit parameter of order k = f(V)
            params.scatter(VFiltered, fit[k], color=colors[k], s=100, alpha=0.7)
            
## Initialize the total fit residual
fitResiduals = 0.

## Loop on the order of the polynomial
for i in range(0, len(paramList[0])):
    
    ## Fit the parameters of the same order for different vacancy numbers
    ## with a polynomial of order VFitOrder
    fit, [res, _, _, _] = poly.polyfit(VList, [row[i] for row in paramList], VFitOrder, full=True)
    ## Get the new fit function
    fitFunc = poly.Polynomial(fit)
    ## Add the residual value to the total one
    fitResiduals = fitResiduals + res
    
    ## Plot the fit and store the parameters
    params.plot(range(1,50), fitFunc(range(1,50)), color=colors[i])
    paramFinalList.append(fit)

## Definition of the obtained 2D fitted function
def resultFunc(x, y):
    f = 0
    for i in range(0, len(paramList[0])):
        for k in range(0, len(paramFinalList[0])):
            f += paramFinalList[i][k] * math.pow(x,i) * math.pow(y,k)
            
    return f

# formationEnergyFit = []
# for i in range(0, len(He)):
#     formationEnergyFit.append(resultFunc(He[i]/V[i], V[i]))
# 
# results.scatter(He/V, V, formationEnergy, c='g', alpha=0.5)
# results.scatter(He/V, V, formationEnergyFit, c='m', alpha=0.5)

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
params.set_title("Fit parameters", fontsize=20)
params.set_xlabel("Vacancy number",fontsize=16)
params.set_xlim([0, 50])
energies.set_title("Formation Energies", fontsize=20)
energies.set_xlabel("Helium/Vacancy number",fontsize=16)

## Print the total residual of the fits
print "Residual of the total fit: ", fitResiduals

## Show the plots
plt.show()