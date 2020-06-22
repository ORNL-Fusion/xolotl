/* =====================================================================================
					 The UQ Toolkit (UQTk) version 3.0
					 Copyright (2015) Sandia Corporation
					 http://www.sandia.gov/UQToolkit/

	 Copyright (2015) Sandia Corporation. Under the terms of Contract
DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains certain
rights in this software.

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
=====================================================================================
*/
/*! \file custom_inf.cpp
 */

#include "arrayio.h"
#include "arraytools.h"
#include "func.h"
#include "inference.h"
#include "mcmc.h"
#include "mrv.h"
#include "post.h"
#include "tools.h"
#include <getopt.h>
#include <unistd.h>

#include <iostream>
#include <map>
#include <sstream>
#include <string>

#include "cmath"

using namespace std;

#define PI 3.14159265359

Array2D<double>
forwardFunc(Array2D<double>& p, Array2D<double>& x, void* funcinfo);

///  Main program: Bayesian inference of a few standard function types
int
main(int argc, char* argv[])
{
	// Input arguments
	int power = 2;
	void* funcinfo = (void*)&power;
	string liktype = "gausmarg";
	string priortype = "uniform";
	double priora = -DBL_MAX;
	double priorb = DBL_MAX;
	Array2D<double> xdata, ydata;
	const char* xfile = "xdata.txt";
	const char* yfile = "Normalized_ydata.txt";
	read_datafileVS(xdata, xfile);
	read_datafileVS(ydata, yfile);
	Array2D<double> xgrid = xdata;
	int nxgr = xgrid.XSize();
	int dataNoiseInference = 1;
	Array1D<double> datanoise(xdata.XSize(), 0.01);
	int pdim = 2;
	int order = 1;
	Array1D<int> rndInd(pdim, 0);
	rndInd(1) = 1;
	string pdftype = "pct";
	int nmcmc = 400000;
	double mcmcgamma = 0.01;
	bool optimflag = true;
	int chdim = 6;
	Array1D<double> chstart(chdim, 2.0);
	Array1D<double> chsig(chdim, 1.0);
	double likParam = 0.000001;
	double likParam_int = 0;
	Array2D<double> pgrid, pchain;
	int nburn = 100000;
	int nstep = 100;

	// Output containers
	Array1D<double> mapparam, pmean_map, pvar_map, fmean_map, fvar_map;
	Array1D<double> datavar_map;
	Array1D<double> p_postave_mean(pdim), p_postave_var(pdim),
		p_postvar_mean(pdim);
	Array1D<double> f_postave_mean(nxgr), f_postave_var(nxgr),
		f_postvar_mean(nxgr);
	Array1D<double> postave_datavar;
	Array2D<double> pmeans, pvars, fmeans, fvars, datavars, paramPCcfs;

	// Run the inference
	infer_model(forwardFunc, funcinfo, liktype, priortype, priora, priorb,
		xdata, ydata, xgrid, dataNoiseInference, datanoise, pdim, order, rndInd,
		pdftype, nmcmc, mcmcgamma, optimflag, chstart, chsig, likParam,
		likParam_int, pgrid, pchain, nburn, nstep, mapparam, datavar_map,
		pmean_map, pvar_map, fmean_map, fvar_map, postave_datavar,
		p_postave_mean, p_postave_var, p_postvar_mean, f_postave_mean,
		f_postave_var, f_postvar_mean, paramPCcfs);

	// Write the outputs

	write_datafile(pchain, "pchain.dat");
	write_datafile_1d(mapparam, "mapparam.dat");

	array1Dto2D(p_postave_mean, pmeans);
	pmeans.insertCol(pmean_map, 1);
	write_datafile(pmeans, "pmeans.dat");
	array1Dto2D(p_postave_var, pvars);
	pvars.insertCol(p_postvar_mean, 1);
	pvars.insertCol(pvar_map, 2);
	write_datafile(pvars, "pvars.dat");

	array1Dto2D(f_postave_mean, fmeans);
	fmeans.insertCol(fmean_map, 1);
	write_datafile(fmeans, "fmeans.dat");
	array1Dto2D(f_postave_var, fvars);
	fvars.insertCol(f_postvar_mean, 1);
	fvars.insertCol(fvar_map, 2);
	write_datafile(fvars, "fvars.dat");

	array1Dto2D(postave_datavar, datavars);
	datavars.insertCol(datavar_map, 1);
	write_datafile(datavars, "datavars.dat");

	write_datafile(paramPCcfs, "parampccfs.dat");

	return 0;
}

Array2D<double>
forwardFunc(Array2D<double>& p, Array2D<double>& x, void* funcinfo)
{
	int* power = (int*)funcinfo;

	int np = p.XSize();
	int nx = x.XSize();
	int pdim = p.YSize();
	int xdim = x.YSize();

	assert(pdim == 2);
	assert(xdim == 1);

	Array2D<double> y(np, nx);

	//  x(0,0)=0.0;
	//   y(0,0)=0.0;

	for (int ip = 0; ip < np; ip++)
		for (int ix = 0; ix < nx; ix++)
			y(ip, ix) = 1 / (p(ip, 0) * x(ix, 0) * sqrt(2 * PI)) *
				exp(-pow(log(x(ix, 0)) - p(ip, 1), *power) /
					(2.0 * p(ip, 0) * p(ip, 0)));

	return y;
}
