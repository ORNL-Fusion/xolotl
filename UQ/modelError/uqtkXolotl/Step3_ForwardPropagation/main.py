import numpy as np
import os
import matplotlib.pyplot

 
# # Append mapparam to pchain
# cmd='cp ../Step1_ParametricFit/Gumbel/pchain.dat .'
# os.system(cmd)
# cmd='cp ../Step1_ParametricFit/Gumbel/mapparam.dat .'
# os.system(cmd)

cmd = "./transpose_file.x mapparam.dat >> pchain.dat"
os.system(cmd)


# Run model inference code in a non-inference, propagation-only regime 
# Parameters have to be very carefully selected to match them with weibull_inference code
METHOD="gausmarg"
DATANOISE=0.0
cmd = 'echo "0" > rndInd_file.dat'
os.system(cmd)
cmd = 'echo "1" >> rndInd_file.dat'
os.system(cmd)


##########################################################################################

std_ModErr, std_PusFor, std_ModErrMAP = np.loadtxt('pvars.txt', usecols = (0,1,2), unpack=True)
p, p_MAP = np.loadtxt('pmeans.txt', usecols = (0,1,), unpack=True)

std_p0 = np.sqrt(std_ModErr[0]+std_PusFor[0])
std_p1 = np.sqrt(std_ModErr[1]+std_PusFor[1])

amid = p[0]
bmid = p[1]
adel = 3.0*std_p0
bdel = 3*std_p1

print "amid = ", amid
print "adel = ", adel
print "bmid = ", bmid
print "bdel = ", bdel

cmd="awk '{print ($1-amid)/adel,$2/adel,($3-bmid)/bdel,$4/bdel,$5/bdel}' amid=" + str(amid) + " adel=" + str(adel) + " bmid=" + str(bmid) + " bdel=" + str(bdel) +  " pchain.dat > pchain_sc.dat"
os.system(cmd)

###########################################################################################

# cmd="awk '{print ($1-amid)/adel,$2/adel,($3-bmid)/bdel,$4/bdel,$5/bdel}' amid=1.0 adel=1.0 bmid=1.0 bdel=0.4 pchain.dat > pchain_sc.dat"
# os.system(cmd)

cmd = "$UQBINDIR/model_inf -z -x Depth.dat -y Depth.dat -f pc -p pchain_sc.dat -l" + str(METHOD) + " -e " + str(DATANOISE) + " -m0 -o 3 -s pct -d 2 -t Depth.dat -r rndInd_file.dat"
os.system(cmd)

# Postprocess
cmd = "paste Depth.dat fmeans.dat fvars.dat > full_propagation.dat"
os.system(cmd)

cmd = "cp full_propagation.dat ../Step4_Postprocessing/full_propagation.dat"
os.system(cmd)

# Clean the rest
cmd = "rm -rf *.dat *.txt"
os.system(cmd)