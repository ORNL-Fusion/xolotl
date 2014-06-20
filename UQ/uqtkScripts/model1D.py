#!/usr/bin/env python

# This file contains the 1D surrogate model for the formation
# energy fit that (forward) propagates the uncertainties.
#
# Use the following to run this script:
#    ./model.py <inputFile> <outputFile>
#
# Ex.  ./model.py quadpts_mapped.dat qd_out.dat
#
# NOTE:  This step is done after pce_quad has been called
#
# TODO 2D

import os
import getopt
import shutil
import sys
import numpy as np
import math
import scipy.special
import scipy.optimize

#############################################################

def model(modelPar):
    return model_example(modelPar)

#############################################################

def model_example(modelPar):

# Evaluates the PC surrogate

    # Load the (scaled to [-1,1]) inputdata
    xin = np.loadtxt('inputdata.dat', usecols = (0,), unpack=True)
    nout=len(xin)
    npar=modelPar.shape[0]
    output=np.empty((npar,nout))
    for j in range(nout):
        order=5        #order of the fit TODO make this a commandline arg
        leg=np.empty((order+1,))
        for iord in range(order+1):
        # Evaluate the Legendre polynomial of order iord at xin[j]
            leg[iord]=scipy.special.eval_legendre(iord,xin[j])
        for i in range(npar):
            print i+1, ": Running the model with parameter setting ", modelPar[i,:]
            # Evaluate the model by taking the dot product of the previously
            # evaluated Legendre polynomial with the PC coeff
            output[i,j]=np.dot(modelPar[i,:],leg)

    return output


#############################################################

def main(argv):
    modelPar_file=argv[0]    # quadpts_mapped.dat
    output_file=argv[1]        # qd_out.dat
    modelPar=np.loadtxt(modelPar_file) # TODO what if 2d?

    output=model(modelPar)
    np.savetxt(output_file,output)

if __name__ == "__main__":
    main(sys.argv[1:])

