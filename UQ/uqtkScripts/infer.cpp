/* =====================================================================================
 The UQ Toolkit (UQTk) version 2.0
 Copyright (2013) Sandia Corporation
 http://www.sandia.gov/UQToolkit/

 Copyright (2013) Sandia Corporation. Under the terms of Contract
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

#include "PCSet.h"
#include "XMLUtils.h"
#include "XMLreader.h"
#include "posterior.h"
#include "uqtkmcmc.h"
#include "uqtktools.h"
#include <math.h>

using namespace std;

/// Main program: MCMC inference of the model parameters given data
int
main(int argc, char* argv[])
{
	// Pointer to posterior information
	postAux* pinfo = new postAux;
	// Read the xml tree
	RefPtr<XMLElement> xmlTree = readXMLTree("infer.xml");
	// Read the model specification and send them to the posterior information
	// structure
	readXMLModelInput(xmlTree, pinfo->modelparams, pinfo->modelparamnames,
		pinfo->modelauxparams);
	// Read specific information needed by inference, e.g. data and the
	// parameters to be inferred
	readXMLDataInput(
		xmlTree, pinfo->data, pinfo->postparams, &(pinfo->noisetype));

	int order = pinfo->modelauxparams(0);
	int dim = pinfo->modelauxparams(1);
	cout << order << " " << dim << endl;

	// Set the PC to use Legendre polynomials
	PCSet pcm("NISPnoq", order, dim, "LU");
	pinfo->pcmodel = &pcm;

	// Number of MCMC steps
	int nsteps;
	// Array to hold the starting values of the chain
	Array1D<double> chstart;
	// Define the MCMC object
	MCMC mchain(LogPosterior, (void*)pinfo);
	// Read the xml file for MCMC-specific information
	readXMLChainInput(xmlTree, &mchain, chstart, &nsteps, pinfo->chainParamInd,
		pinfo->priortype, pinfo->priorparam1, pinfo->priorparam2);

	// Prepend the parameter names to the output file
	FILE* f_out;
	string filename = mchain.getFilename();
	int chdim = chstart.XSize();

	f_out = fopen(filename.c_str(), "w");

	fprintf(f_out, "%s ", "Step");
	for (int i = 0; i < chdim; i++)
		fprintf(f_out, "%21s ",
			pinfo->modelparamnames(pinfo->chainParamInd(i)).c_str());
	fprintf(f_out, "%24s %24s \n", "Accept_prob", "Log_posterior");
	fclose(f_out);

	// Set the flag
	mchain.namesPrepended();
	// Run the chain
	mchain.runOptim(chstart);
	mchain.runChain(nsteps, chstart);

	return 0;
}
