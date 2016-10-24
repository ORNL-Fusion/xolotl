import math
import numpy as np
import matplotlib.pyplot as plt
import os

x = np.loadtxt('xdata.txt', usecols = (0,), unpack=True)
y = np.loadtxt('Normalized_ydata.txt', usecols = (0,), unpack=True)
mu_mean, mu_MAP = np.loadtxt('fmeans.dat', usecols = (0,1), unpack=True)
sigma_m, sigma_p, sigma_m_MAP = np.loadtxt('fvars.dat', usecols = (0,1,2), unpack=True)
sigma_inf = np.loadtxt('mapparam.dat', usecols = (0,), unpack=True)

# # Uncomment for combined data
# x_plot = np.zeros((len(mu_mean)/2))
# mu_mean_plot = np.zeros((len(mu_mean)/2))
# sigma_m_plot = np.zeros((len(mu_mean)/2))
# sigma_p_plot = np.zeros((len(mu_mean)/2))
# sigma_m_MAP_plot = np.zeros((len(mu_mean)/2))
# for i in range(len(x_plot)):
#     x_plot[i] = x[i]
#     mu_mean_plot[i] = mu_mean[i]
#     sigma_m_plot[i] = sigma_m[i]
#     sigma_p_plot[i] = sigma_p[i]
#     sigma_m_MAP_plot[i] = sigma_m_MAP[i]

#Uncomment for single type data (MD or SRIM)
x_plot = x
mu_mean_plot = mu_mean
sigma_m_plot = sigma_m
sigma_p_plot = sigma_p
sigma_m_MAP_plot = sigma_m_MAP

# # Uncomment if posterior predictive is inferred as log(sigma)
# sigma_PostPred = math.exp(sigma_inf[5])

# # Uncomment if posterior predictive is inferred as sigma
# sigma_PostPred = sigma_inf[5]

# Uncomment if posterior predictive is not inferred
sigma_PostPred = 0

PFP1 = np.sqrt(sigma_p_plot)
PFP2 = np.sqrt(sigma_p_plot+sigma_m_plot)
PFP3 = np.sqrt(sigma_p_plot+sigma_m_plot+sigma_PostPred**2)

PFP_MAP1=np.sqrt(sigma_m_MAP_plot)
PFP_MAP2=np.sqrt(sigma_m_MAP_plot+sigma_PostPred**2)

fig = plt.figure(figsize=(10,7))
ax=fig.add_axes([0.10,0.15,0.85,0.75])
plt.fill_between(x_plot,mu_mean_plot-PFP3,mu_mean_plot+PFP3,color='cyan',label='Posterior predictive')
plt.fill_between(x_plot,mu_mean_plot-PFP2,mu_mean_plot+PFP2,color='blue',label='Model error')
plt.fill_between(x_plot,mu_mean_plot-PFP1,mu_mean_plot+PFP1,color='black',label='Pushed forward posterior')

# #Uncomment to plot MAP values
# plt.fill_between(x_plot,mu_mean_plot-PFP_MAP2,mu_mean_plot+PFP_MAP2,color='cyan',label='Post predictive stdev')
# plt.fill_between(x_plot,mu_mean_plot-PFP_MAP1,mu_mean_plot+PFP_MAP1,color='blue',label='Pushed forward stdev + model error')

plt.plot(x_plot, mu_mean_plot, color='red', linewidth=2, label='Mean prediction')
plt.plot(x,y,'o', color='red',label='Data')
ax.set_xlabel("Depth (nm)",fontsize=22)
ax.set_ylabel("PDF",fontsize=22)
# plt.ylim((0,0.25))
plt.legend(loc='best')
plt.savefig('Gumbel_UQ.jpg')
plt.clf()

# # Save the important results
cmd = "cp pvars.dat ../../Step2_SurrogateConstruction/pvars.txt"
os.system(cmd)
cmd = "cp pmeans.dat ../../Step2_SurrogateConstruction/pmeans.txt"
os.system(cmd)
cmd = "cp pchain.dat ../../Step3_ForwardPropagation/pchain.dat"
os.system(cmd)
cmd = "cp mapparam.dat ../../Step3_ForwardPropagation/mapparam.dat"
os.system(cmd)
cmd = "cp pvars.dat ../../Step3_ForwardPropagation/pvars.txt"
os.system(cmd)
cmd = "cp pmeans.dat ../../Step3_ForwardPropagation/pmeans.txt"
os.system(cmd)
