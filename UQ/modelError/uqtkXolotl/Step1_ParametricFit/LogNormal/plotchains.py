import numpy as np
import matplotlib.pyplot as plt

# #Uncomment if posterior predictive is also inferred
# chain0, chain1, chain2, chain3, chain4, chain5 = np.loadtxt('pchain.dat', usecols = (0,1,2,3,4,5), unpack=True)

#Uncomment if posterior predictive is not inferred
chain0, chain1, chain2, chain3, chain4 = np.loadtxt('pchain.dat', usecols = (0,1,2,3,4), unpack=True)

MCMC_step = np.zeros((len(chain0)))

for i in range(len(chain0)):
    MCMC_step[i]=i

fig = plt.figure(figsize=(10,7))
ax=fig.add_axes([0.10,0.15,0.85,0.75])
plt.plot(MCMC_step, chain0, color='red', linewidth=2, label='Mean prediction')
ax.set_xlabel("MCMC Step",fontsize=22)
ax.set_ylabel("Parameter 0",fontsize=22)
plt.legend(loc='best')
plt.savefig('Chain0.jpg')
plt.clf()

fig = plt.figure(figsize=(10,7))
ax=fig.add_axes([0.10,0.15,0.85,0.75])
plt.plot(MCMC_step, chain1, color='red', linewidth=2, label='Mean prediction')
ax.set_xlabel("MCMC Step",fontsize=22)
ax.set_ylabel("Parameter 1",fontsize=22)
plt.legend(loc='best')
plt.savefig('Chain1.jpg')
plt.clf()

fig = plt.figure(figsize=(10,7))
ax=fig.add_axes([0.10,0.15,0.85,0.75])
plt.plot(MCMC_step, chain2, color='red', linewidth=2, label='Mean prediction')
ax.set_xlabel("MCMC Step",fontsize=22)
ax.set_ylabel("Parameter 2",fontsize=22)
plt.legend(loc='best')
plt.savefig('Chain2.jpg')
plt.clf()

fig = plt.figure(figsize=(10,7))
ax=fig.add_axes([0.10,0.15,0.85,0.75])
plt.plot(MCMC_step, chain3, color='red', linewidth=2, label='Mean prediction')
ax.set_xlabel("MCMC Step",fontsize=22)
ax.set_ylabel("Parameter 3",fontsize=22)
plt.legend(loc='best')
plt.savefig('Chain3.jpg')
plt.clf()

fig = plt.figure(figsize=(10,7))
ax=fig.add_axes([0.10,0.15,0.85,0.75])
plt.plot(MCMC_step, chain4, color='red', linewidth=2, label='Mean prediction')
ax.set_xlabel("MCMC Step",fontsize=22)
ax.set_ylabel("Parameter 4",fontsize=22)
plt.legend(loc='best')
plt.savefig('Chain4.jpg')
plt.clf()

# #Uncomment if posterior predictive is also inferred
# fig = plt.figure(figsize=(10,7))
# ax=fig.add_axes([0.10,0.15,0.85,0.75])
# plt.plot(MCMC_step, chain5, color='red', linewidth=2, label='Mean prediction')
# ax.set_xlabel("MCMC Step",fontsize=22)
# ax.set_ylabel("STD",fontsize=22)
# plt.legend(loc='best')
# plt.savefig('ChainSTD.jpg')
# plt.clf()