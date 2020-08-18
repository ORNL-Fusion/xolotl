#!/usr/bin/env python
import numpy as np
import math
import matplotlib.pyplot as plt
from   pylab import *


## Create plots
fig = plt.figure()
ax3 = plt.subplot(111)

## Open the data file from Xolotl (by default retentionOut.txt)
fluence1, He1, D1, T1, V1, I1 = loadtxt('/home/sophie/Data/Xolotl/LANL/retention_4e25.txt', usecols = (0,1,2,3,4,5) , unpack=True)

## Helium retention
ax3.plot(fluence1, 100.0 * He1/fluence1, lw=6, ls='-', color="black", alpha=0.5, label='He')

ax3.set_xlabel("Fluence [nm-2]",fontsize=25)
ax3.set_ylabel("Retention [%]",fontsize=25)
#ax3.set_xlim([0.0, 10.5])
#ax3.set_ylim([0.0, 0.145])
#ax3.set_yscale('log')
ax3.tick_params(axis='both', which='major', labelsize=25)
ax3.tick_params(axis='both', which='minor', labelsize=25)
## Plot the legend
l = ax3.legend(loc='best')
setp(l.get_texts(), fontsize=25)

## Show the plots
plt.show()
