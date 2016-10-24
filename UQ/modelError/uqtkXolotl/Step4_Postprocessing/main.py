import math
import numpy as np
import matplotlib.pyplot as plt
import os

depth, cum_mean, cum_MAP, sigma_m, sigma_p, sigma_m_MAP  = np.loadtxt('full_propagation.dat', usecols = (0,1,2,3,4,5), unpack=True)  

PFP1 = np.zeros(len(depth))
PFP2 = np.zeros(len(depth))

for i in range(len(depth)):
    a = max(sigma_p[i]+sigma_m[i],0.0)
    PFP2[i] = math.sqrt(a)
    b =  max(sigma_p[i],0.0)
    PFP1[i] = math.sqrt(b)
#     sigma_m_MAP[i] = math.sqrt(sigma_m_MAP[i])
# 
# PFP1 = np.sqrt(sigma_p)
# PFP2 = np.sqrt(sigma_p+sigma_m)

fig = plt.figure(figsize=(12,7))
ax=fig.add_axes([0.10,0.15,0.85,0.75])
plt.fill_between(depth,cum_mean-PFP2,cum_mean+PFP2,color='blue',label='Model error')
plt.fill_between(depth,cum_mean-PFP1,cum_mean+PFP1,color='black',label='Pushed forward posterior')
plt.plot(depth, cum_mean, color='red', linewidth=2, label='Mean prediction')
ax.set_xlabel("Depth (nm)",fontsize=22)
# ax.set_ylabel("Cumulative He1 Fraction",fontsize=22)
# ax.set_ylabel("Cumulative He1 Distribution",fontsize=22)
ax.set_ylabel("Cumulative Fraction of Retained Helium",fontsize=22)
plt.ylim((0.0,1.0))
plt.legend(loc=4)
plt.savefig('UQ_result.jpg')
plt.clf()

# Clean the rest
cmd = "rm -rf *.dat *.txt"
os.system(cmd)