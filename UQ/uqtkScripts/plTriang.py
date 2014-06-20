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
#
# Given a set of samples of random variables, this script plots a lower triangular
# matrix of marginalized densities. The diagonal contains the density of individual
# random variables, marginalized over all other variables. Plots below the diagonal
# contain the 2D density of the associated pair of random variables, marginalized over
# all other variables.
#
# The script expects an input file containing samples of random variable samples. Each
# column contains the samples of a different random variable. Each line is a sample
# of the full random variable set.
#
# The first line of the file needs to contain labels for the respective columns
#
# The script is written for handling Markov Chain Monte Carlo (MCMC) posterior sample files,
# so it assumes that the first column is an integer MCMC step number, and that the last two columns
# contain the acceptance probability (labeled "acceptance_prob") and posterior probability. The script
# can easily be modified to process files without those columns.
#
# Note 1) If the samples come from a full MCMC file, then the column labeled
#   acceptance_prob will have a value of -1 on the last line. This line
#   contains the MAP values of the chain, and must be removed from
#   the sample list as it is not a regular sample. If the samples
#   come from a partial MCMC run (e.g. a run in progress), then
#   all lines need to be kept.
# Note 2) All columns should be of equal length
#
# Written spring 2013 by Cosmin Safta and Bert Debusschere, Sandia National Labs, Livermore, CA
# <bjdebus@sandia.gov>

import os
import sys
import getopt
import math

import numpy as npy
import matplotlib.pyplot as plt
from scipy import stats, mgrid, c_, reshape, random, rot90
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
from scipy.stats import norm

#sys.path.append('../../math_libs/pytools')
import file_utils

help_string = """
Usage:
  plot_triang.py [-h] -i <mcmc filename> [--nb <burn-in samples>]  [-s <stride>] [-d] [--nkde <n_pts_kde>]
what
  Plot MCMC variable marginal and covariance PDF graphs
where
  -h = print help info
  -i = name of file containing MCMC data
  -s = stride with which to read the file, default 1
  --nb = number of burn-in samples to be removed from the chain, default 0
  -d = use dense format: plot full set of plots without any labels, white space. Print out
       set of figures for each pair of 2 variables only (useful when more than 5 or 6 MCMC variables).
       default = false
  --nkde = number of points per dimension to evaluate densities at with KDE, default = 100
           (more points gives a smoother PDF but increases the computational cost)
""" 
# Set to 1 to get more output to screen
# Set to > 1 to get a lot of output to screen
debug = 0

# Set to 1 for showing plots interactively
interact = 0


#
# Process inputs 
#
try:
  opts, extra  = getopt.getopt(sys.argv[1:],"hi:s:d",["nb=","nkde="])
except getopt.GetoptError, err:
    print str(err)
    print help_string
    sys.exit(1)

if (debug > 0):
    print "Extra command line arguments following recognized arguments: ",extra

# Default values
samples_file_name=""
n_burnin = 0
stride = 1
np_kde=100

for o,a in opts:
  if o == "-h":
    print help_string
    sys.exit(0)
  elif o == "-i":
    samples_file_name = a
  elif o == "-s":
    stride = int(a)
  elif o == "--nb":
    n_burnin = int(a)
  elif o == "--nkde":
    np_kde = int(a)
  else:
    assert False, "unhandled command line parsing option"

# error checking
if(samples_file_name==""):
  print "Sample file name must be specified"
  print help_string
  sys.exit(1)

if (n_burnin < 0):
  print "The number of burn-in samples needs to be >= 0"
  print help_string
  sys.exit(1)

if (stride < 0): 
  print "The file read stride needs to be >= 0"
  print help_string
  sys.exit(1)

if (np_kde < 1): 
  print "The number of KDE points per dimension needs to be >= 1"
  print help_string
  sys.exit(1)

#
# Import data from MCMC file
#
# Section 1
# load chain file
# it expects first column is line id, then the last
# two columns are alpha and current log posterior 
print "Loading in data file",samples_file_name
all_samples, vnames = file_utils.extract_all_vars(samples_file_name,n_burnin,debug,stride)
n_all_vars = len(vnames)
n_cols = len(all_samples[0,:])

# Extract all MCMC chain variables in separate array
d0 = all_samples[:,1:1+n_all_vars]
if (debug > 0):
    print d0.shape

# Some settings to connect with code Cosmin gave me
nthin = 1                # take only every nthin state (for faster kde)
nskip = 0                # entries to skip
#npdf = 100               # no of grid points for 1D pdf's
istart = 0               # number of column with first MCMC variable
#cend   = 0               # number of extra columns at the end.
cend   = 0               # remove columns at end so code runs faster for debugging

nvars=d0.shape[1]-istart-cend  # number of variables we will actually process
print 'Number of sample lines in file',d0.shape[0]
print 'Number of vars  we will process in file',nvars

# Section 2
# set up 2D kde objects
print "Setting up 2D KDE objects"
kern_i_j=[]
for j in range(istart+1,istart+nvars):
    for i in range(istart,j):
        if (debug > 2):
            print i,j
        kern_i_j.append(stats.kde.gaussian_kde(c_[d0[nskip::nthin,i],d0[nskip::nthin,j]].T))

# Section 3
# set up 2D meshes and evaluate kde objects on those meshes
# no. of grid points is controlled with kde_idx, defaults to 100
print "Evaluating 2D KDE objects on meshes. This may take a while ..."
kde_idx = np_kde*1j # complex number to include end points
xmesh=[]; ymesh=[]; zmesh=[];
icount=0
cov_idx = npy.zeros((nvars,nvars),dtype=npy.int) # 2D array to keep track of which index in xmesh etc. the 
                                                   # the plots corresponding to vars i,j belong to
for j in range(istart+1,istart+nvars):
    for i in range(istart,j):
        if (debug > 0):
            print "Computing 2D marginal distribution between variables:",i,",",j,":",vnames[i]," & ",vnames[j]
        x,y = mgrid[d0[nskip:,i].min():d0[nskip:,i].max():kde_idx, d0[nskip:,j].min():d0[nskip:,j].max():kde_idx]
        z   = reshape(kern_i_j[icount](c_[x.ravel(), y.ravel()].T).T, x.T.shape)
        xmesh.append(x);
        ymesh.append(y);
        zmesh.append(z);
        cov_idx[i,j] = icount
        icount = icount+1

# Section 4
# evaluate 1D pdfs
print "Evaluating 1D marginal pdfs with KDE and Gaussian fit"
xlin=[]; pdflin=[]; pdfgauss=[];
for i in range(istart,istart+nvars):
    xlin.append(npy.linspace(d0[nskip:,i].min(),d0[nskip:,i].max(),np_kde)) ;
    kernlin=stats.kde.gaussian_kde(d0[nskip::nthin,i]);
    mean, deviation = norm.fit(d0[nskip::nthin,i]);
    pdflin.append(kernlin(xlin[i-istart]));
    pdfgauss.append(norm.pdf(xlin[i-istart],mean,deviation));
    
    print i, mean, deviation



# Section 5
print "Assembling tri-diagonal plots in non-dense format"

# ds is the distance between subplots
# xs,ys are the coordinates (normalized) of the subplot in the lower left corner
# xe,ye are the distances left in the uppper right corner
# fsizex, fsizey are figure sizes
# ncont are no of contours for 2D pdfs
xs=0.12; ys=0.1; ds=0.04
xe=0.08; ye=0.05
fsizex=12; fsizey=12;
ncont=20; 
sx=(1-(nvars-1)*ds-xs-xe)/nvars;
sy=(1-(nvars-1)*ds-ys-ye)/nvars;
fs1=20
majorFormatter = FormatStrFormatter('%6.0e')

fig = plt.figure(figsize=(fsizex,fsizey))

# Section 5.1
subs=[]
# add diagonal plots
for i in range(nvars):
    subs.append(fig.add_axes([xs+i*(sx+ds),ys+(nvars-1-i)*(sy+ds),sx,sy]))

# add lower triangular plots
for i in range(nvars-1):
    for j in range(i+1):
        if (debug > 2):
            print j,(nvars-2-i)
        subs.append(fig.add_axes([xs+j*(sx+ds),ys+(nvars-2-i)*(sy+ds),sx,sy]))

subsnp=npy.array(subs)

# Plot 1D pdfs
for i in range(nvars):
    subsnp[i].plot(xlin[i],pdflin[i])
    # Gaussians
    subsnp[i].plot(xlin[i],pdfgauss[i],c=(1,0,0))      

# Plot 2D pdfs
for i in range(nvars*(nvars-1)/2):
    subsnp[nvars+i].contour(xmesh[i],ymesh[i],zmesh[i],ncont)

# Section 5.2
# no ticks and ticklabels
for subpl in subsnp:
    subpl.set_xticks([]); 
    subpl.set_yticks([]); 

# for diagonal plots only put xmin and xmax
for i in range(istart,istart+nvars):
    subsnp[i-istart].set_xticks([d0[nskip:,i].min(),d0[nskip:,i].max()]); 

# Set y labels on the right for diagonal plots
for i in range(nvars):
    subsnp[i].yaxis.tick_right()
    subsnp[i].yaxis.set_label_position("right")
    subsnp[i].set_ylabel(vnames[i], fontsize=fs1)

    plt.savefig("posteriors.pdf")