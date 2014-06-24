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

# Input settings 
# They have to be changed to the value they have in the XML file
xorder=5
yorder=3

chdimlist=[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25]      # Indices of chain parameters
chstart=[69.734004, 50.788312, -6.280163, 0.674329, 49.214234, 37.250268, -3.675781, 0.671578, 20.428427, 18.915853, -0.516389, 0.329204, -1.647776, -0.460337, 0.766950, -0.208465, -1.370898, -1.764382, -0.413855, 0.006777, 0.507299, 0.650158, 0.241819, 0.108867, -0.737367]  # Starting value of the chain
chainfile='chain.dat'      # Name of the chain file
n_burnin=0              # Burnin samples

## Prepare the xml file ##################################################
shutil.copyfile('infer2D.xml.templ','infer2D.xml')
for line in fileinput.input('infer2D.xml', inplace = 1):
    print line.replace('CHAINFILE', chainfile),



xin, yin, zin = loadtxt('inputdata2d.dat', usecols = (0,1,2), unpack=True)

xin=2.*(xin/3.8125)-1.
yin=2.*((yin-1)/43.)-1.
np.savetxt("inputdata.dat", np.transpose([xin,yin,zin]))

## Run the inference code #####################################################
ii=len(chdimlist)
for idim in chdimlist:
    ii=ii-1
    for line in fileinput.input("infer2D.xml", inplace = 1):
        print line.replace('PAR_'+str(ii+1), str(chstart[ii])),
 
print "Running the parameter inference"
os.system('./infer2D.x')


## Import data from MCMC file ###################################################
print "Loading in chain file",chainfile
all_samples, vnames = file_utils.extract_all_vars(chainfile,n_burnin,0,1)
n_all_vars = len(vnames)
n_cols = len(all_samples[0,:])
# Extract all MCMC chain variables in separate array
chn = all_samples[:,0:1+n_all_vars]
nchn=chn.shape[0]


## Plot chains ##################################################################
for i in range(n_all_vars):
    fig = plt.figure(figsize=(10,7))
    ax=fig.add_axes([0.10,0.15,0.85,0.75])

    plt.plot(chn[:,0],chn[:,i+1],color='black',linewidth=2)
    ax.set_xlabel("MCMC step",fontsize=22)
    ax.set_ylabel(vnames[i],fontsize=22)

    plt.savefig('chn_'+vnames[i]+'.pdf')
    plt.clf()

for i in range(n_all_vars):
    for j in range(i):
        fig = plt.figure(figsize=(10,7))
        ax=fig.add_axes([0.10,0.15,0.85,0.75])
    
        plt.plot(chn[:,j+1],chn[:,i+1],'ko',markeredgecolor='black',markersize=5)
        ax.set_xlabel(vnames[j],fontsize=22)
        ax.set_ylabel(vnames[i],fontsize=22)

        plt.savefig('chn_'+vnames[j]+'_'+vnames[i]+'.pdf')
        plt.clf()


## Find posterior predictive #######################################################
# ideally we should have a deterministic forward model code,
# but for a line it is simple enough to compute here
print "Sampling chain to compute posterior predictive distribution"
param = np.zeros(len(chdimlist))

ngr=len(yin)
zpp=np.zeros((ngr))
z2psh=np.zeros((ngr))
zstd=np.zeros((ngr))
for ip in range(nchn):
    ii=0
    for idim in chdimlist:
        # grab variabl samples from posterior chain
        param[idim-1]=chn[ip,ii+1]
        ii=ii+1

    # run forward model for these sampled parameter values
    
    for j in range(ngr):
        leg = np.empty(len(chdimlist)-1)
        i=0
        for a in range(xorder+1):
            for b in range(yorder+1):
                leg[i]=scipy.special.eval_legendre(a,xin[j])*scipy.special.eval_legendre(b,yin[j])
                i=i+1

        newparam=param[:-1]
        zval=np.dot(newparam,leg)
        # running sum for y
        zpp[j]=zpp[j]+zval
        # running sum for y^2
        z2psh[j]=z2psh[j]+zval**2
        # running sigma
        zstd[j]=zstd[j]+exp(param[len(chdimlist)-1])

# pushed forward mean
zpp=zpp/nchn
# pushed forward ave y^2
z2psh=z2psh/nchn

zstd=zstd/nchn
# pushed forward std dev
std_push=(z2psh-zpp**2)**0.5
# Posterior predictive std dev
std_pp=(std_push**2 + zstd**2)**0.5# Scale back x and y
xin = (xin+1.) * 3.8125/2.
yin = ((yin+1.) * 43./2.) + 1.
# and convert y to integers
yin = np.asarray([int(y+0.5) for y in yin])

# Loop on the list of V numbers
Vlist = {1, 2, 6, 14, 18, 19, 27, 32, 44}

# Plot pushed forward posterior and posterior predictive 
fig = plt.figure(figsize=(10,7))
ax=fig.add_axes([0.10,0.15,0.85,0.75])
for v in Vlist:
  # Filter everything
  xinF = xin[yin == v]
  zinF = zin[yin == v]
  zppF = zpp[yin == v]
  std_ppF = std_pp[yin == v]
  std_pushF = std_push[yin == v]
  yinF = yin[yin == v]
  
  plt.fill_between(xinF,zppF-std_ppF,zppF+std_ppF,color='lightgrey',label='Post predictive stdev')
  plt.fill_between(xinF,zppF-std_pushF,zppF+std_pushF,color='grey',label='Pushed forward stdev')
  plt.plot(xinF, zppF, linewidth=1, color='red', label='Mean prediction')
  plt.plot(xinF, zinF,'o', markersize=4, color='black', label='Data')

ax.set_xlabel("x",fontsize=22)
ax.set_ylabel("y",fontsize=22)

plt.savefig('postpred.pdf')
plt.clf()

# Plot pushed forward posterior and posterior predictive difference
fig = plt.figure(figsize=(10,7))
ax=fig.add_axes([0.10,0.15,0.85,0.75])
for v in Vlist:
  # Filter everything
  xinF = xin[yin == v]
  zinF = zin[yin == v]
  zppF = zpp[yin == v]
  std_ppF = std_pp[yin == v]
  std_pushF = std_push[yin == v]
  yinF = yin[yin == v]
  
  plt.fill_between(xinF,-std_ppF,std_ppF,color='lightgrey',label='Post predictive stdev')
  plt.fill_between(xinF,-std_pushF,std_pushF,color='grey',label='Pushed forward stdev')
  plt.plot(xinF, zppF-zppF, linewidth=1, color='red', label='Mean prediction')
  plt.plot(xinF, zinF-zppF,'o', markersize=4, color='black', label='Data')

ax.set_xlabel("x",fontsize=22)
ax.set_ylabel("y",fontsize=22)

plt.savefig('postpreddiff.pdf')
plt.clf()


## Plot posterior 'triangle' ####################################################
# calling it as a script for now
cmd='./plTriang2D.py -i ' + chainfile + ' --nb ' + str(n_burnin) + ' > plTriang.log'
print "Running",cmd
os.system(cmd)

print "END: bayesianInference2D done: check out *.pdf files for results."
