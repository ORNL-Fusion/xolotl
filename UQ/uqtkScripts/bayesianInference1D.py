#!/usr/bin/env python
#=====================================================================================
#                     The UQ Toolkit (UQTk) version 2.0
#                    Copyright (2013) Sandia Corporation
#                     http://www.sandia.gov/UQToolkit/
#
#    Copyright (2013) Sandia Corporation. Under the terms of Contract DE-AC04-94AL85000
#    with Sandia Corporation, the U.S. Government retains certain rights in this software.
#
#    This file is part of The UQ Toolkit (UQTk)
#
#    UQTk is free software: you can redistribute it and/or modify
#    it under the terms of the GNU Lesser General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    UQTk is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU Lesser General Public License for more details.
#
#    You should have received a copy of the GNU Lesser General Public License
#    along with UQTk.  If not, see <http://www.gnu.org/licenses/>.
#
#    Questions? Contact Bert Debusschere <bjdebus@sandia.gov>
#    Sandia National Laboratories, Livermore, CA, USA
#====================================================================================


import os
import shutil
import sys
import numpy as np
import math
import random as rnd
from scipy.stats.mstats import mquantiles
import scipy.optimize
import matplotlib.pyplot as plt

import fileinput
import file_utils

from pylab import *

rc('legend',loc='upper right', fontsize=20)
rc('lines', linewidth=4, color='r')
rc('axes',linewidth=3,grid=True,labelsize=22)
rc('xtick',labelsize=16)
rc('ytick',labelsize=16)

###########################################################################
# define uqtkbin
if os.environ.get("UQTK_SRC") is None:
    print "Error: Need to set path to uqtk src as environment variable UQTK_SRC -> Abort"
    quit()
else:
    if ( not os.path.isdir(os.environ["UQTK_SRC"]) ):
        print "\"",os.environ["UQTK_SRC"],"\" is not a valid path -> Abort"
        quit()

uqtkbin=os.environ["UQTK_SRC"]+"/src_cpp/bin"
pcequad=uqtkbin+"/pce_quad"
###########################################################################

## Input settings 

# Indices of chain parameters
chdimlist=[1,2,3,4,5,6,7] 
# Starting value of the chain (or of the optimization if it is used)    
chstart=np.ones((7))  
# Name of the chain file                 
chainfile='chain.dat'  
# Burn-in samples     
n_burnin=3000          

## Prepare the xml file ##################################################
shutil.copyfile('infer.xml.templ','infer.xml')
for line in fileinput.input('infer.xml', inplace = 1):
    print line.replace('CHAINFILE', chainfile),

## Scale the x between -1 and 1 ################################################
# Open the original file
xin, yin = loadtxt('inputdata1d.dat', usecols = (0,1), unpack=True)
# Scale (here for V = 27)
xin=2.*(xin/4.44444)-1.
# Save in the file that will be used for the MCMC
np.savetxt("inputdata.dat", np.transpose([xin,yin]))

## Run the inference code #####################################################
ii=0
for idim in chdimlist:
    for line in fileinput.input("infer.xml", inplace = 1):
        print line.replace('PAR_'+str(idim), str(chstart[ii])),
    ii=ii+1

print "Running the parameter inference"
os.system('./infer.x')

## Import data from MCMC file ###################################################
print "Loading in chain file",chainfile
all_samples, vnames = file_utils.extract_all_vars(chainfile,n_burnin,0,1)
n_all_vars = len(vnames)
n_cols = len(all_samples[0,:])
# Extract all MCMC chain variables in separate array
chn = all_samples[:,0:1+n_all_vars]
nchn=chn.shape[0]

## Find posterior predictive #######################################################
# ideally we should have a deterministic forward model code,
# but for a line it is simple enough to compute here
print "Sampling chain to compute posterior predictive distribution"
param = np.zeros(len(chdimlist))

ngr=len(yin)
ypp=np.zeros((ngr))
y2psh=np.zeros((ngr))
ystd=np.zeros((ngr))
for ip in range(nchn):
    ii=0
    for idim in chdimlist:
        # grab variabl samples from posterior chain
        param[idim-1]=chn[ip,ii+1]
        ii=ii+1

    # run forward model for these sampled parameter values
    xlin=np.linspace(xin[0],xin[len(xin)-1],ngr)
    
    for j in range(ngr):
        order=len(chdimlist)-2# or -1
        leg=np.empty((order+1,))
        for iord in range(order+1):
            leg[iord]=scipy.special.eval_legendre(iord,xlin[j])

        yval=np.dot(param[:-1],leg)
        # running sum for y
        ypp[j]=ypp[j]+yval
        # running sum for y^2
        y2psh[j]=y2psh[j]+yval**2
        # running sigma
        ystd[j]=ystd[j]+exp(param[order+1])

# pushed forward mean
ypp=ypp/nchn
# pushed forward ave y^2
y2psh=y2psh/nchn
# inferred standard deviation
ystd=ystd/nchn
# pushed forward std dev
std_push=(y2psh-ypp**2)**0.5
# Posterior predictive std dev
std_pp=(std_push**2 + ystd**2)**0.5

# Scale back the x
xlin = (xlin + 1.) * 4.44444/2.

# Plot pushed forward posterior and posterior predictive 
fig = plt.figure(figsize=(10,7))
ax=fig.add_axes([0.10,0.15,0.85,0.75])
plt.fill_between(xlin,ypp-std_pp,ypp+std_pp,color='lightgrey',label='Post predictive stdev')
plt.fill_between(xlin,ypp-std_push,ypp+std_push,color='grey',label='Pushed forward stdev')
plt.plot(xlin, ypp, linewidth=1, label='Mean prediction')
plt.plot(xlin, yin,'o', markersize=4, color='black',label='Data')
ax.set_xlabel("x",fontsize=22)
ax.set_ylabel("y",fontsize=22)
plt.legend()

plt.savefig('postpred.pdf')
plt.clf()

# Plot pushed forward posterior and posterior predictive difference
fig = plt.figure(figsize=(10,7))
ax=fig.add_axes([0.10,0.15,0.85,0.75])
plt.fill_between(xlin,ypp-std_pp-ypp,ypp+std_pp-ypp,color='lightgrey',label='Post predictive stdev')
plt.fill_between(xlin,ypp-std_push-ypp,ypp+std_push-ypp,color='grey',label='Pushed forward stdev')
plt.plot(xlin, ypp-ypp, linewidth=1, label='Mean prediction')
plt.plot(xlin, yin-ypp,'o', markersize=4, color='black',label='Data')
ax.set_xlabel("x",fontsize=22)
ax.set_ylabel("y",fontsize=22)
plt.legend()

plt.savefig('postpreddiff.pdf')
plt.clf()


## Plot chains ##################################################################
# 1D
# Comment down to plt.clf() if you don't want the 1D chains
for i in range(n_all_vars):
    fig = plt.figure(figsize=(10,7))
    ax=fig.add_axes([0.10,0.15,0.85,0.75])

    plt.plot(chn[:,0],chn[:,i+1],color='black',linewidth=2)
    ax.set_xlabel("MCMC step",fontsize=22)
    ax.set_ylabel(vnames[i],fontsize=22)

    plt.savefig('chn_'+vnames[i]+'.pdf')
    plt.clf()

# 2D
# Comment down to plt.clf() if you don't want the 2D chains
for i in range(n_all_vars):
    for j in range(i):
        fig = plt.figure(figsize=(10,7))
        ax=fig.add_axes([0.10,0.15,0.85,0.75])
      
        plt.plot(chn[:,j+1],chn[:,i+1],'ko',markeredgecolor='black',markersize=5)
        ax.set_xlabel(vnames[j],fontsize=22)
        ax.set_ylabel(vnames[i],fontsize=22)
  
        plt.savefig('chn_'+vnames[j]+'_'+vnames[i]+'.pdf')
        plt.clf()


## Plot posterior 'triangle' ####################################################
# calling it as a script for now
# Comment down to os.system(cmd) if you don't want the triangle plot
cmd='./plTriang.py -i ' + chainfile + ' --nb ' + str(n_burnin) + ' > plTriang.log'
print "Running",cmd
os.system(cmd)

print "END: Bayesian Inference done: check out *.pdf files for results."
