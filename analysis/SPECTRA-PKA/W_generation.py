## Author: Jose March-Rico
## Revisions: Nick Militello, Sophie Blondel

import numpy as np
from   pylab import *

## Constants
N_A = 6.02214129e23
density = 19.254 #g cm-3
atomic_weight = 183.84
max_size = 200

## Read the PKA data 
lo, hi, pka = loadtxt('W_fast_elemental.dat', usecols = (0,1,2) , unpack=True)
pka_E = (pka * 0.001 * N_A * density)/atomic_weight # PKA energy in keV
pka_rate_cubic = (lo + hi) * 1.0e-21 * 5.0e5 # PKA/s/nm^3

## For the following tables, the columns give fractions of clusters of size n produced
filenamev = "W_vac.txt"
filenamei = "W_int.txt"

## Read in the vacancy production table:
casc_E = [0] # PKA energy for cascade
defect_frac_v = [] # Fractions for defect production. n x m array (n based on PKA energy, m based on defect size)
defect_row = [0] * max_size
defect_frac_v.append(defect_row) # To consider E = 0 keV case (no production)
with open(filenamev,"r") as f:
    for line in f:
        array = list(map(float,line.strip().split(" ")))
        temp = array[1:]
        temp.extend([0.0] * (max_size - len(temp)))
        defect_frac_v.append(temp) # fractions for each size
        casc_E.append(float(array[0])) # energy

## Read in the interstitial production table (energies are the same as for V):
defect_frac_i = [] # Fractions for defect production. n x m array (n based on PKA energy, m based on defect size)
defect_frac_i.append(defect_row) #To consider E = 0 keV case (no production)
with open(filenamei,"r") as f:
    for line in f:
        array = list(map(float,line.strip().split(" ")))
        temp = array[1:]
        temp.extend([0.0] * (max_size - len(temp)))
        defect_frac_i.append(temp) # fractions for each size

## The defect flux calculation will be done for each pka_E in the spectrum, and based on the cubic PKA rate (s^-1 nm^-3)
defect_flux_cubic_v = [0] * max_size
defect_flux_cubic_i = [0] * max_size
avgFracI = []
avgFracV = []
ei1 = 0
ei2 = 0

## For each energy in the spectrum, add the expected defect flux to the total:
for pe in range(1,len(pka_E)):

    ## Find the energy indices in casc_E to perform interpolation (note that pka_E is always > 0)
    for ce in range(ei1,len(casc_E)):
        ## If the PKA energy is greater than the greatest data point (40 keV)
        if(pka_E[pe] > casc_E[-1]): 
            ## Create an averaged defect fraction vector based on sub-cascade 1.0 scaling:
            avgFracI = []
            avgFracV = []

            ## The rows are for PKA energies 
            ## Columns are defect sizes
            for dd in range(0, max_size): # The index doesn't matter here
                # Average the last five energy bins
                tempI, tempV = 0.0, 0.0
                for i in range(1,6):
                    tempI += defect_frac_i[-i][dd]*(pka_E[pe]/casc_E[-i])
                    tempV += defect_frac_v[-i][dd]*(pka_E[pe]/casc_E[-i])
                    
                avgFracI.append(tempI / 5.0)
                avgFracV.append(tempV / 5.0)
    
        ## Get the indexes for energy pe
        elif(casc_E[ce] > pka_E[pe]):
            ei1 = ce-1 # first index for energy pe
            ei2 = ce # second index for energy pe
            break


    ## Now add the defect flux from cascades of the given PKA energy using linear interpolation between data points:
    ## Integrate along the PKA energy spectrum:
    for di in range(0,len(defect_frac_i[0])):

        #Interstitials:
        #Find the defect fluxes at energy pe and pe-1, using linear interpolation:
        if (pka_E[pe] <= casc_E[-1]): 
            defect_flux_1 = (defect_frac_i[ei1][di] + ((defect_frac_i[ei2][di]-defect_frac_i[ei1][di])/(casc_E[ei2]-casc_E[ei1]))*(pka_E[pe]-casc_E[ei1]))*pka_rate_cubic[pe]
        else: 
            defect_flux_1 = avgFracI[di]*pka_rate_cubic[pe]
        if (defect_flux_1 > 0.0): 
            defect_flux_cubic_i[di] += defect_flux_1
        
        if (pka_E[pe] <= casc_E[-1]): 
            defect_flux_1 = (defect_frac_v[ei1][di] + ((defect_frac_v[ei2][di]-defect_frac_v[ei1][di])/(casc_E[ei2]-casc_E[ei1])*(pka_E[pe]-casc_E[ei1])))*pka_rate_cubic[pe]
        elif (pka_E[pe] > casc_E[-1]): defect_flux_1 = avgFracV[di]*pka_rate_cubic[pe]
        if(defect_flux_1 > 0.0): 
            defect_flux_cubic_v[di] += defect_flux_1
        
#Make sure that the rate of interstitial production = rate of vacancy production:
sumI = 0
sumV = 0
for ii in range(0,len(defect_frac_i[0])):
    sumI += defect_flux_cubic_i[ii]*(ii+1)
    sumV += defect_flux_cubic_v[ii]*(ii+1)
sumI *= 1.0e21
sumV *= 1.0e21
    
print("Comparing the total rate of interstitials and vacancies:")
print("{:9e}".format(sumI))
print("{:9e}".format(sumV))
print()

factor = sumV / sumI
#Increase the interstitial production rate until it equals the vacancy rate:
for ff in range(0,len(defect_flux_cubic_i)): 
    defect_flux_cubic_i[ff] = defect_flux_cubic_i[ff] * factor
sumI = 0
for ii in range(0,len(defect_frac_i[0])):
    sumI+= defect_flux_cubic_i[ii]*(ii+1)
sumI *= 1.0e21
print("Corrected rate of interstitials and vacancies:")
print("{:e}".format(sumI))
print("{:e}".format(sumV))

#Print final results:
print()
print("Interstitial fluxes:")
print(defect_flux_cubic_i)
print()
print("Vacancy fluxes:")
print(defect_flux_cubic_v)


## Create plots
fig1 = plt.figure()
pkaPlot = plt.subplot(111)

## Create the x axis
x = []
for i in range(max_size):
    x.append(i)
    
    
pkaPlot.scatter(x, defect_flux_cubic_i, s=50, color='magenta', alpha=1, label='Interstitial')
pkaPlot.scatter(x, defect_flux_cubic_v, s=50, color='black', alpha=1, label='Vacancy')

## Plot the legend
l = pkaPlot.legend(loc='best')
setp(l.get_texts(), fontsize=20)

## Some shaping
pkaPlot.set_xlabel("Number of Defects per Cluster",fontsize=25)
pkaPlot.set_ylabel("Rate [nm$^{-3}$ s$^{-1}$]",fontsize=25)
pkaPlot.set_xlim([0, 160])
#pkaPlot.set_ylim([1.0e6, 1.0e13])
#pkaPlot.set_xscale('log')
pkaPlot.set_yscale('log')
pkaPlot.tick_params(axis='both', which='major', labelsize=25)
pkaPlot.tick_params(axis='both', which='minor', labelsize=25)

## Show the plots
plt.show()
