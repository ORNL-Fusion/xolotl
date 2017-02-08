/* =====================================================================================
                     The UQ Toolkit (UQTk) version 3.0
                     Copyright (2015) Sandia Corporation
                     http://www.sandia.gov/UQToolkit/

     Copyright (2015) Sandia Corporation. Under the terms of Contract DE-AC04-94AL85000
     with Sandia Corporation, the U.S. Government retains certain rights in this software.

     This file is part of The UQ Toolkit (UQTk)

     UQTk is free software: you can redistribute it and/or modify
     it under the terms of the GNU Lesser General Public License as published by
     the Free Software Foundation, either version 3 of the License, or
     (at your option) any later version.

     UQTk is distributed in the hope that it will be useful,
     but WITHOUT ANY WARRANTY; without even the implied warranty of
     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
     GNU Lesser General Public License for more details.

     You should have received a copy of the GNU Lesser General Public License
     along with UQTk.  If not, see <http://www.gnu.org/licenses/>.

     Questions? Contact Bert Debusschere <bjdebus@sandia.gov>
     Sandia National Laboratories, Livermore, CA, USA
===================================================================================== */
/*! \file custom_inf.cpp
*/

#include <unistd.h>
#include <sstream>
#include <map>
#include <iostream>
#include <string>

#include <getopt.h>

#include "cmath"

#include "func.h"
#include "post.h"
#include "mrv.h"
#include "inference.h"

#include "mcmc.h"
#include "tools.h"
#include "arrayio.h"
#include "arraytools.h"



using namespace std;

Array2D<double> forwardFunc(Array2D<double>& p, Array2D<double>& x, void* funcinfo);



//  Main program: Bayesian inference of model parameters
int main (int argc, char *argv[]) 
{
    // Input arguments
    // Any auxiliary information for model can be passed via void*
    int power=2;
    void* funcinfo=(void*) &power;
    // Likelihood type
    string liktype="gausmarg";
    // Prior type  
    string priortype="uniform";
    // Prior bounds
    double priora=0.0;
    double priorb=DBL_MAX;
    // Read data
    Array2D<double> xdata,ydata;
    const char* xfile="xdata.txt";
    const char* yfile="Normalized_ydata.txt";
    read_datafileVS(xdata,xfile);
    read_datafileVS(ydata,yfile);
   // Prediction points
    Array2D<double> xgrid=xdata; 
    int nxgr=xgrid.XSize();

    // Indicates whether data noise is inferred or not
    int dataNoiseInference=1;
    // Data noise standard deviation, either fixed, or initial position 	(if inferred)
    Array1D<double> datanoise(xdata.XSize(),0.01);
    // Number of parameters of the model
    int pdim=2;
    // Order for prediction PC
    int order=1;
    // Indices of parameters to be 'randomized'
    Array1D<int> rndInd(pdim,0); rndInd(1)=1;
    // PDF type of parameters
    string pdftype="pct";
    // PC type of parameters
    string pctype="HG";
    // MCMC seed
    int seed=13;
    // MCMC starting point
    int nmcmc=500000;
    // Gamma factor for AMCMC
    double mcmcgamma=0.0001;
    // Flag whether to prepend optimization or not
    bool optimflag=true;
    // Chain dimensionality
    int chdim=6;
    // Number of burnin steps
    int nburn=100000;
    // For good statistics, pick every nstep state
    int nstep=100;
    // Chain start 
    Array1D<double> chstart(chdim,1.0);
    // Standard deviation per dimension for initial non-adaptive part
    Array1D<double> chsig(chdim,0.1);
    // Likelihood parameter (double)
    double likParam=0.000001;
    // Likelihood parameter (int)
    double likParam_int=0;
    // Arrays for parameter values, and clean chain
    Array2D<double> pgrid,pchain;
 
    // Output containers
    Array1D<double> mapparam,pmean_map,pvar_map, fmean_map,fvar_map;
    Array1D<double> datavar_map;
    Array1D<double> p_postave_mean(pdim), p_postave_var(pdim), p_postvar_mean(pdim);
    Array1D<double> f_postave_mean(nxgr), f_postave_var(nxgr), f_postvar_mean(nxgr);
    Array1D<double> postave_datavar;
    Array2D<double> pmeans,pvars,fmeans,fvars,datavars,paramPCcfs;

    
    // Run the inference
    infer_model(forwardFunc,funcinfo,
        liktype,
        priortype,priora,priorb,
        xdata,ydata, xgrid,
        dataNoiseInference,datanoise,
        pdim,order,rndInd,pdftype,pctype,
        seed, nmcmc,mcmcgamma,optimflag,chstart,chsig,
        likParam,likParam_int,
        pgrid, pchain, nburn, nstep,
        mapparam, datavar_map,
        pmean_map,pvar_map, 
        fmean_map,fvar_map,
        postave_datavar,
        p_postave_mean,p_postave_var,p_postvar_mean, 
        f_postave_mean,f_postave_var,f_postvar_mean,
        paramPCcfs);

    // Write the outputs
    
      write_datafile(pchain,"pchain.dat");
      write_datafile_1d(mapparam,"mapparam.dat");    
      
      array1Dto2D(p_postave_mean, pmeans);
      pmeans.insertCol(pmean_map,1);
      write_datafile(pmeans,"pmeans.dat");
      array1Dto2D(p_postave_var, pvars);
      pvars.insertCol(p_postvar_mean,1);
      pvars.insertCol(pvar_map,2);
      write_datafile(pvars,"pvars.dat");

      array1Dto2D(f_postave_mean, fmeans);
      fmeans.insertCol(fmean_map,1);
      write_datafile(fmeans,"fmeans.dat");
      array1Dto2D(f_postave_var, fvars);
      fvars.insertCol(f_postvar_mean,1);
      fvars.insertCol(fvar_map,2);
      write_datafile(fvars,"fvars.dat");

      array1Dto2D(postave_datavar,datavars);    
      datavars.insertCol(datavar_map,1);
      write_datafile(datavars,"datavars.dat");

      write_datafile(paramPCcfs,"parampccfs.dat");
   
  return 0;
}

Array2D<double> forwardFunc(Array2D<double>& p, Array2D<double>& x, void *funcinfo)
{
 //   int* power=(int*) funcinfo;

    int np=p.XSize();
    int nx=x.XSize();
    int pdim=p.YSize();
    int xdim=x.YSize();

    assert(pdim==2);
    assert(xdim==1);
    
    Array2D<double> y(np,nx,0.0);

    for(int ip=0;ip<np;ip++)
	{
//	printf("\n%lg %lg\n",p(ip,0),p(ip,1));
	if (p(ip,0)<=0.0)
		return y;
	if (p(ip,1)<=0.0)
		return y;

        for(int ix=0;ix<nx;ix++)
            y(ip,ix)=p(ip,0)/p(ip,1) * pow(x(ix,0)/p(ip,1),(p(ip,0)-1.0)) * exp(-pow(x(ix,0)/p(ip,1),p(ip,0)));
	}
    return y;
}
