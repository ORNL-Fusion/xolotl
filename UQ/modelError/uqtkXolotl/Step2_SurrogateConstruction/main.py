import numpy as np
import os
import matplotlib.pyplot

import cPickle as pick

# Build the domain file for Weibull/Gumbel parameters
# Need to pick a and b range around the region where the MCMC is concentrated
# Select the domain margins for both k and lambda
std_ModErr, std_PusFor, std_ModErrMAP = np.loadtxt('pvars.txt', usecols = (0,1,2), unpack=True)
p, p_MAP = np.loadtxt('pmeans.txt', usecols = (0,1,), unpack=True)

std_p0 = np.sqrt(std_ModErr[0]+std_PusFor[0])
std_p1 = np.sqrt(std_ModErr[1]+std_PusFor[1])

lambdaMargin0 = p[0] - 3*std_p0
lambdaMargin1 = p[0] + 3*std_p0
kMargin0 = p[1] - 3*std_p1
kMargin1 = p[1] + 3*std_p1

lambdaMargin0=max(0.0,lambdaMargin0)
kMargin0=max(0.0,kMargin0)

# Construct the domain
cmd = 'echo ' + str(lambdaMargin0) + " " + str(lambdaMargin1) + ' > pdomain.dat' # lambda range
os.system(cmd)
 
cmd = 'echo ' + str(kMargin0) + " " + str(kMargin1) + ' >> pdomain.dat' # k range
os.system(cmd)
 
# Select the order of the surrogate
ORDER = 3
# Typically set points-per-dimension that is at least one higher than the order
PPD = 5
 
# Prepare input file (training points) for model runs
# See uq_pc.py, or run 'uq_pc.py -h' for details 
cmd = "$UQPCDIR/uq_pc.py -r offline_prep -p pdomain.dat -m proj -s quad -n "+str(PPD)+" -v 33"
os.system(cmd)
 
# Run Xolotl at training points
p0, p1 = np.loadtxt('ptrain.dat', usecols = (0,1), unpack=True)
  
# for i in range(1):
for i in range(len(p0)):
    cmd = "echo " + str(p0[i]) + " " + str(p1[i]) + " > FitParameters.dat"
    os.system(cmd) 
    cmd = "mpiexec -n 4 /home/ocekmer/Workspaces/Sourceforge/uq-build/xolotl /home/ocekmer/Workspaces/Sourceforge/uq-build/params.txt"
    os.system(cmd)

#     W_Depth, CumulHe = np.loadtxt('heliumCumul_92.dat', usecols = (0,1), unpack=True)
    Fluence, Reten = np.loadtxt('retentionOut.txt', usecols = (0,1,), skiprows=1, unpack=True)
#     CumulHe = CumulHe / max(CumulHe)
#     np.savetxt('TrainOut.dat',CumulHe)
    np.savetxt('TrainOut.dat',Reten)
    if i==0:
        cmd = "./transpose_file.x TrainOut.dat > ytrain.dat"
        os.system(cmd)
    else:
        cmd = "./transpose_file.x TrainOut.dat >> ytrain.dat"
        os.system(cmd)

# np.savetxt('Depth.dat',W_Depth)
np.savetxt('Depth.dat',Fluence)

# Run Xolotl at validation points
p0val, p1val = np.loadtxt('pval.dat', usecols = (0,1), unpack=True)
 
for i in range(len(p0val)):
    cmd = "echo " + str(p0val[i]) + " " + str(p1val[i]) + " > FitParameters.dat"
    os.system(cmd) 
    cmd = "mpiexec -n 4 /home/ocekmer/Workspaces/Sourceforge/uq-build/xolotl /home/ocekmer/Workspaces/Sourceforge/uq-build/params.txt"
    os.system(cmd)

#     W_DepthVal, CumulHeVal = np.loadtxt('heliumCumul_92.dat', usecols = (0,1), unpack=True)
    FluenceVal, RetenVal = np.loadtxt('retentionOut.txt', usecols = (0,1,), skiprows=1, unpack=True)
#     CumulHeVal = CumulHeVal / max(CumulHeVal)
#     np.savetxt('TrainOut.dat',CumulHeVal)
    np.savetxt('TrainOut.dat',RetenVal)
    if i==0:
        cmd = "./transpose_file.x TrainOut.dat > yval.dat"
        os.system(cmd)
    else:
        cmd = "./transpose_file.x TrainOut.dat >> yval.dat"
        os.system(cmd)
     
# Compute the surrogate of order $ORDER
cmd = "$UQPCDIR/uq_pc.py -r offline_post -p pdomain.dat -m proj -s quad -n " + str(PPD) + " -v 33 -t " + str(ORDER)
os.system(cmd)
 
# Plot data-versus-model for surrogate accuracy assessment
cmd = "$UQPCDIR/plot.py dm training validation"
os.system(cmd)
 
results=pick.load(open('results.pk', 'rb'))

# ndgrid=len(W_Depth)
ndgrid=len(Fluence)

mi=np.loadtxt('mindex.dat')
npc=mi.shape[0]
pccf_all=np.empty((npc,ndgrid))
for i in range(ndgrid):
    pccf_all[:,i]=results['pcmi'][0][i]
 
np.savetxt('pccf_all.dat',pccf_all)
 
# Save the important results
cmd = "cp Depth.dat ../Step3_ForwardPropagation/Depth.dat"
os.system(cmd)
 
cmd = "cp pccf_all.dat ../Step3_ForwardPropagation/"
os.system(cmd)
 
cmd = "cp mindex.dat ../Step3_ForwardPropagation/mindexp.dat"
os.system(cmd)
 
# Clean the rest
cmd = "rm -rf *.dat *.log *.png pcf *.pk"
os.system(cmd)