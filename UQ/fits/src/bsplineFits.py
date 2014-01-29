#!/usr/bin/env python
#======================================================================================
#Created on Jan 28, 2014
# Code to fit binding energy of Helium as a function of Helium and Vacancy cluster sizes
# using Smooth and Least Squares Bivariate Spline (bspline) approximations 
#@author: Crystal Jernigan
#=======================================================================================

import numpy as np
import math
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from   pylab import *
from   scipy import stats
from scipy import interpolate

# Load data from 'HeVexpEHe.dat'
# x = col[0], y = col[1], z = col[2]
x, y, z = loadtxt('HeVexpEHe.dat', usecols = (0,1,2) , unpack=True)   

# Smooth bivariate spline approximation
# kx, ky are the degrees of the bivariate spline in the x and y direction respectively
zsfit5 = interpolate.SmoothBivariateSpline(x,y,z,kx=5,ky=5)
# Spline approximation evaluated at x,y
speval5 = zsfit5.ev(x,y)

# Least-squares bivariate spline approximation with specified 1D sequences of knot coordinates
xknots1 = array([0,1,240,249,250])
yknots1 = array([0,1,40,49,50])
lsqzs1fit5 = interpolate.LSQBivariateSpline(x,y,z,xknots1.flat,yknots1.flat,kx=5,ky=5)
lsqev1 = lsqzs1fit5.ev(x,y)

xknots2 = array([0,1,40,75,183,194,222,231,250])
yknots2 = array([0,1,5,11,43,45,46,47,50])
lsqzs2fit5 = interpolate.LSQBivariateSpline(x,y,z,xknots2.flat,yknots2.flat,kx=5,ky=5)
lsqev2 = lsqzs2fit5.ev(x,y)

xknots3 = array([0,1,14,125,165,184,219,220,223,225,231,232,239,241,250])
yknots3 = array([0,1,11,24,25,30,43,45,46,47,48,49,50])
lsqzs3fit5 = interpolate.LSQBivariateSpline(x,y,z,xknots3.flat,yknots3.flat,kx=5,ky=5)
lsqev3 = lsqzs3fit5.ev(x,y)

#PLOT 3D FIGURE
fig = plt.figure()
rect = fig.add_subplot(2,2,1).get_position()
ax = Axes3D(fig, rect)
ax.set_ylim3d([0,10])
ax.set_xlabel('Helium Size', fontsize=12)
ax.set_ylabel('Vacancy Size', fontsize=12)
ax.set_zlabel('exp(-E_He)', fontsize=12)
ax.set_title('Scaled Binding Energy of Helium vs. Helium, Vacancy Cluster Sizes', fontsize=14)
# Plot of input data
ax.scatter(x,y,z,zdir='z',s=16, color='g')

# Plots of fits to exact data
rect = fig.add_subplot(2,2,2).get_position()
ax = Axes3D(fig, rect)
ax.set_xlabel('Helium Size', fontsize=12)
ax.set_ylabel('Vacancy Size', fontsize=12)
ax.set_zlabel('exp(-E_He)', fontsize=12)
ax.set_title('Smooth Bivariate Spline Fit, kx=5, ky=5', fontsize=14)
ax.scatter(x,y,speval5,zdir='z')
ax.scatter(x,y,z,zdir='z',s=16,color='g')

rect = fig.add_subplot(2,2,3).get_position()
ax = Axes3D(fig, rect)
ax.set_xlabel('Helium Size', fontsize=12)
ax.set_ylabel('Vacancy Size', fontsize=12)
ax.set_zlabel('exp(-E_He)', fontsize=12)
ax.set_title('LSQ Bivariate Spline Fit, kx=5, ky=5; knots1', fontsize=14)
ax.scatter(x,y,lsqev1,zdir='z', color='r')
ax.scatter(x,y,z,zdir='z',s=16,color='g')

rect = fig.add_subplot(2,2,4).get_position()
ax = Axes3D(fig, rect)
ax.set_xlabel('Helium Size', fontsize=12)
ax.set_ylabel('Vacancy Size', fontsize=12)
ax.set_zlabel('exp(-E_He)', fontsize=12)
ax.set_title('LSQ Bivariate Spline Fit, kx=5, ky=5; knots3', fontsize=14)
ax.scatter(x,y,z,zdir='z', color='g')
ax.scatter(x,y,lsqev3,zdir='z',s=16, color='r')
plt.show()