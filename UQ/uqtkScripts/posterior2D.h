/* =====================================================================================
 The UQ Toolkit (UQTk) version 2.0
 Copyright (2013) Sandia Corporation
 http://www.sandia.gov/UQToolkit/

 Copyright (2013) Sandia Corporation. Under the terms of Contract DE-AC04-94AL85000
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

struct postAux {
	Array2D<double> data;
	Array1D<double> modelparams;
	Array1D<string> modelparamnames;
	Array1D<double> modelauxparams;
	Array1D<double> postparams;
	string noisetype;
	Array1D<string> priortype;
	Array1D<double> priorparam1;
	Array1D<double> priorparam2;
	Array1D<int> chainParamInd;
	PCBasis* pcmodel;
};

/// \brief Evaluate the log of the posterior for a given set of model
/// and nuisance parameters
double LogPosterior(Array1D<double>& m, void *info);

double LogPosterior(Array1D<double>& m, void *info) {
	double pi = 4. * atan(1.);

	postAux* pinfo = (postAux*) info;

	double logprior = 0;

	for (int ic = 0; ic < m.Length(); ic++) {

		if (!strcmp(pinfo->priortype(ic).c_str(), "uniform")) {
			double a = pinfo->priorparam1(ic);
			double b = pinfo->priorparam2(ic);
			if (a >= b)
				throw Tantrum("Prior bounds are incorrect!");
			if (m(ic) < b && m(ic) > a)
				logprior -= log(b - a);
			else
				return -1.e180;
		}

		else if (!strcmp(pinfo->priortype(ic).c_str(), "normal")) {
			double mu = pinfo->priorparam1(ic);
			double sig = pinfo->priorparam2(ic);

			logprior -= 0.5 * log(2. * pi * sig * sig);
			logprior -= 0.5 * pow((m(ic) - mu) / sig, 2.0);
		} else
			throw Tantrum("Only unifrom or normal priors are implemented!");
	}

	Array2D<double> data = pinfo->data;
	// Set the parameter values for the forward run
	int chaindim = m.XSize();
	Array1D<double> modelparams(chaindim - 1);
	for (int ic = 0; ic < chaindim - 1; ic++)
		modelparams(pinfo->chainParamInd(ic)) = m(ic);

	// Posterior parameter
	// Either Proportionaility constant between signal and noise for likelihood construction
	// or standard deviation itself
	double stdpar = pinfo->postparams(0);

	// Compare the model data with measurement data according to the noise model
	// to get the likelihood.
	double std;

	double sum = logprior;
	int nTot = data.XSize();
	int nDim = data.YSize() - 1;

	for (int it = 0; it < nTot; it++) {
		Array1D<double> xx(nDim, 0.e0);
		for (int id = 0; id < nDim; id++) {
			xx(id) = data(it, id);
		}

		int xOrder = pinfo->modelauxparams(0);
		int yOrder = pinfo->modelauxparams(1);

		Array1D<double> xBasis(xOrder + 1, 0.e0);
		Array1D<double> yBasis(yOrder + 1, 0.e0);
		pinfo->pcmodel->EvalBasis(xx(0), xBasis);
		pinfo->pcmodel->EvalBasis(xx(1), yBasis);

		int i = 0;
		double modelval = 0.0;
		for (int a = 0; a < xOrder + 1; a++) {
			for (int b = 0; b < yOrder + 1; b++) {
				modelval += modelparams(i) * xBasis(a) * yBasis(b);
				i++;
			}
		}

		if (!strcmp(pinfo->noisetype.c_str(), "const_stn")) {
			std = stdpar * fabs(modelval);
		} else if (!strcmp(pinfo->noisetype.c_str(), "const_stdev")) {
			std = stdpar;
		} else if (!strcmp(pinfo->noisetype.c_str(), "infer_stdev")) {
			std = exp(m(chaindim - 1));
		} else
			throw Tantrum("Noise type is not recognized!");

		double var = std * std;

		// weights for 2D inference
		Array1D<double> weight(50, 1.0);
		weight(1) = 1.0 / 10.0;
		weight(2) = 1.0 / 15.0;
		weight(6) = 1.0 / 32.0;
		weight(14) = 1.0 / 68.0;
		weight(18) = 1.0 / 88.0;
		weight(19) = 1.0 / 91.0;
		weight(27) = 1.0 / 121.0;
		weight(32) = 1.0 / 122.0;
		weight(44) = 1.0 / 190.0;

		double err = data(it, nDim) - modelval;
		err = pow(err, 2) * weight(xx(1));
		sum = sum - 0.5 * log(2 * pi) - 0.5 * log(var) - err / (2.0 * var);

	}

	return sum;

}
