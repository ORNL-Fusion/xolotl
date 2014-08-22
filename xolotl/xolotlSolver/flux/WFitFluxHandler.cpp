#include "WFitFluxHandler.h"
#include <iostream>
#include <cmath>

using namespace xolotlSolver;

double WFitFluxHandler::FitFunction100(double x) {
	// Value at which the flux goes to 0
	double x1 = 10.0;

	if (x > x1) return 0.0;

	// Compute the fit
	double value = 7.00876507 + 0.6052078 * x - 3.01711048 * pow(x, 2)
	+ 1.36595786 * pow(x, 3) - 0.295595 * pow(x, 4) + 0.03597462 * pow(x, 5)
	- 0.0025142 * pow(x, 6) + 0.0000942235 * pow(x, 7) - 0.0000014679 * pow(x, 8);

	return value;
}

void WFitFluxHandler::initializeFluxHandler(int numGridpoints, double step) {

	// Set the step size
	stepSize = step;

	double normFactor = 0.0;
	for (int i = 1; i < numGridpoints; i++) {
		double x = (double) i * stepSize;

		normFactor += FitFunction100(x) * stepSize;
	}

//	std::cout << normFactor << std::endl;

	// Factor the incident flux will be multiplied by
	double heFluxNormalized = heFlux / normFactor;

	// The first value should always be 0.0 because of boundary conditions
	incidentFluxVec.push_back(0.0);

	// Starts a i = 1 because the first value was already put in the vector
	for (int i = 1; i < numGridpoints; i++) {
		auto x = i * stepSize;

		auto incidentFlux = heFluxNormalized * FitFunction100(x);

		incidentFluxVec.push_back(incidentFlux);
	}

}

