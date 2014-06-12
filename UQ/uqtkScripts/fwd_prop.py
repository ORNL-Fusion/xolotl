#!/usr/bin/env python

# This file uses forward UQ to propagate the uncertainties.
# This file also allows for seperate preprocessing and postprocessing


import os
import getopt
import shutil
import sys
import numpy as np
import math
import random as rnd
from scipy import stats, mgrid, reshape, random
import cPickle as pick

from model1D import model

#############################################################
def usage():
    # TODO more informative usage()
    # Usage string
    usage_str='uq_pc.py -r <run_regime> -i <inpc_file> -p <pc_type> -d <in_pcdim> -o <in_pcord> -q <nqd> -s <sp_type> -t <out_pcord>'
    def_str='uq_pc.py -r online -i pccf_all.dat -p HG -d 2 -o 1 -q 7 -s full -t 3'
    print "Correct syntax is"
    print usage_str
    print "Default values are"
    print def_str
    #print "For more detailed information on parameters, refer to the UQTk Manual"

#######################################################################################

def main(argv):
    
    # Defaults
    run_regime="online"  # Running regime, "online", "offline_prep" or "offline_post"
    inpc_file="pccf_all.dat" # Input PC coefficient file
    pc_type="HG"            # PC type
    in_pcdim=2              # Input PC dimensionality
    in_pcord=1              # Input PC order
    nqd=7                      # Number of quadrature points per dim (if sp_type=full), or level (if sp_type=sparse)
    sp_type="full"              # Quadrature sparsity type, "full" or "sparse"
    out_pcord=3                    # Output PC order
    #nval=0                     # Number of validation samples, uniform random
    
    
    # Flags for input checks
    rflag=False
    iflag=False
    pflag=False
    dflag=False
    oflag=False
    qflag=False
    sflag=False
    tflag=False
    
        
    try:
        opts, args = getopt.getopt(argv,"hr:i:p:d:o:q:s:t:",["regime=","pcfile=","pctype=","pcdim=","pcord=","nqd=","sparsity=","outord"])
    except getopt.GetoptError:
        usage()
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            usage()
            sys.exit()
        elif opt in ("-r", "--regime"):
            run_regime = arg
            rflag=True
        elif opt in ("-i", "--pcfile"):
            inpc_file = arg
            iflag=True
        elif opt in ("-p", "--pctype"):
            pc_type = arg
            pflag=True
        elif opt in ("-d", "--pcdim"):
            in_pcdim = int(arg)
            dflag=True
        elif opt in ("-o", "--pcord"):
            in_pcord = int(arg)
            oflag=True
        elif opt in ("-q", "--nqd"):
            nqd = int(arg)
            qflag=True
        elif opt in ("-s", "--sparsity"):
            sp_type = arg
            sflag=True
        elif opt in ("-t", "--outord"):
            out_pcord = int(arg)
            tflag=True
      
    # Load parameter domain file
    if (os.path.isfile(inpc_file)):
        pcf_all=np.atleast_2d(np.loadtxt(inpc_file))
        print pcf_all
    else:
        print "Error: The requested input PC coefficient file %s does not exist. Exiting." % inpc_file
        sys.exit()

    # Get the dimensions
    npar=pcf_all.shape[1]
    #npc=pcf_all.shape[0] # not used

    # TODO Sanity checks, e.g. on npc
    #print "Error: The domain file %s contains wrong bounds. Check the row number %d. Exiting." % (pdomain_file,i+1)
    #sys.exit()

    # Print the inputs for reference
    print "Run regime                        ", run_regime
    print "Input PC coefficient file       ", inpc_file
    print "PC type                           ", pc_type
    print "Input PC dim                      ", in_pcdim
    print "Input PC order                    ", in_pcord
    print "Sparsity type                     ", sp_type
    print " with parameter                   ", nqd
    print "Output PC order                   ", out_pcord

    # (1) Generate sample points for online or offline_prep regimes 
    if run_regime=="online" or run_regime=="offline_prep":
        # TODO LU or CC?
    # TODO is this part needed?  since we are already doing the pce_quad
        cmd="generate_quad -d"+str(in_pcdim)+"  -g"+pc_type+" -x"+sp_type+" -p"+str(nqd)+" > gq.log"
        print "Running "+cmd
        os.system(cmd)

        inqdp=np.loadtxt('qdpts.dat')
        npt=inqdp.shape[0]
        print "Quadrature points are in %s in a format %d x %d " % ('qdpts.dat',npt,in_pcdim)

        # Evaluate input PCs at quadrature points
        np.savetxt('xdata.dat',inqdp)
        inpar=np.empty((npt,npar))
        for i in range(npar):
            np.savetxt('pcf',pcf_all[:,i])
            cmd="pce_eval -x PC -s"+pc_type+" -o"+str(in_pcord)+" -f 'pcf' > pcev.log"
            print "Running "+cmd
            os.system(cmd)
            inpar[:,i]=np.loadtxt('ydata.dat')
        
        np.savetxt('qd_in.dat',inpar)
        # Cleanup
        cmd="rm -rf indices.dat"
        os.system(cmd)

        # Exit if only sample preparation is required
        if run_regime=="offline_prep":
            print "Preparation of samples is done."
            sys.exit()
                
    ############################################################################
                
    # (2) Load sample points for online or offline_post regimes
    inpar=np.loadtxt('qd_in.dat').reshape(-1,npar)
    #if (nval>0):
    #    pval_unsc=np.loadtxt(input_val).reshape(-1,dim)

    npt=inpar.shape[0]
    print "Number of quadrature points for forward propagation : "+str(npt)

    ############################################################################

    # (3) Get model outputs

    # Run the model online or....
    if run_regime=="online":
        output=model(inpar)
        #if (nval>0):
        #    yval=model(pval_unsc)

    # ...or read the results from offline simulations
    elif run_regime=="offline_post":
        output=np.atleast_2d(np.loadtxt('qd_out.dat'))
        #if (nval>0):
        #    yval=np.atleast_2d(np.loadtxt("yval.dat"))

    # Read the number of output observables or the number of values of deisgn parameters (e.g. location, time etc..) 
    nout=output.shape[1]
    # Force this for now 
    #assert(nout==1)
    print "Number of output observables of the model is ", nout

    ############################################################################

    # (4) Obtain the PC surrogate using model simulations

    # Emplty arrays and lists to store results
    pccf_all=[]
    mindex_all=[]

    output_pc=np.empty((npt,nout))
    allsens=np.empty((in_pcdim,nout)) #in_pcdim or npar?
    
    evalPCsurrFile = open('evalPCsurr.dat', 'w')
    xin, yin = np.loadtxt('inputdata.dat', usecols = (0,1), unpack=True)
    # Loop over all output observables/locations
    for i in range(nout):
        
        ################################
        
        # (4a) Build PC surrogate
        print "##################################################"
        print "Building surrogate for observable %d / %d" % (i+1,nout)
        np.savetxt('ydata.dat',output[:,i])

        # NOTE: only LEG is implemented in pce_resp
        cmd="pce_resp -x"+pc_type+" -o"+str(out_pcord)+" -d"+str(in_pcdim)+" -e > pcr.log"
        print "Running "+cmd
        os.system(cmd)

        # Get the PC evaluations
        output_pc[:,i]=np.loadtxt('ydata_pc.dat')
        # Get the PC coefficients and multiindex
        pccf=np.loadtxt('PCcoeff_quad.dat')
        mindex=np.loadtxt('mindex.dat')
        

        # Append the results
        pccf_all.append(pccf)
        mindex_all.append(mindex)
    np.savetxt('pccf_'+str(i)+'.dat',pccf)

        ################################
        
        # (4c) Compute moments
        cmd="pce_sens -m'mindex.dat' -f'PCcoeff_quad.dat' -x"+pc_type+" > pcsens.log"
        #cmd="pce_sens -m'mindex.dat' -f'PCcoeff_quad.dat' -x"+pc_type+" > pcsens_"+str(i)+".log"
        print "Running "+cmd
        os.system(cmd)
        
        cmd='grep -e "Mean, Var = " pcsens.log | cut -c 13- > mean_var.log'
        os.system(cmd)
        ins = open("mean_var.log", "r")
        line = ins.readline().split()
        line[0] = line[0] +"\n"
        evalPCsurrFile.write("%s %s %s" %(xin[i],line[0]))
        allsens[:,i]=np.loadtxt('mainsens.dat')
        print "AAA", allsens.shape

    evalPCsurrFile.close()
        
    ############################################################################

    results = {'feval':(inpar,output,output_pc),'pcmi':(pccf_all,mindex_all),'sens':(allsens)}

    # Save results
    pick.dump(results,open('results.pk','wb'),-1)


    # Cleanup of unneeded leftovers
    #del_cmd='rm -rf ydata_pc.dat ydata.dat xdata.dat pccf.dat varfrac.dat totsens.dat mainsens.dat jointsens.dat sp_mindex*dat mindex.dat PCcoeff.dat xwghts.dat wghts.dat qdpts.dat'
    #os.system(del_cmd)


if __name__ == "__main__":
   main(sys.argv[1:])