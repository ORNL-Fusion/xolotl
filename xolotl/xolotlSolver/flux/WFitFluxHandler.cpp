#include "WFitFluxHandler.h"
#include <iostream>
#include <cmath>

using namespace xolotlSolver;

void WFitFluxHandler::initializeFluxHandler(int numGridpoints, double step) {

	// Set the step size
	stepSize = step;

	double x1 = (29.0 - sqrt(41)) / 4.0;
	double normFactor = 1.5 * pow(x1, 4) - 29.0 * pow(x1, 3) + 150 * pow(x1, 2);

	for (int i = 0; i < numGridpoints; i++) {
		auto x = i * stepSize;

		// NOTE: Integrating between the zeros, x0 = 0 and x1 = (29 - sqrt(41))/4,
		// gives the normalizing factor:
		// 1.5 * x^4 - 29 * x^3 + 150 * x^2, where x = (29 - sqrt(41))/4
		auto incidentFlux = (heFlux / normFactor) * (6.0 * x * x * x - 87.0 * x * x + 300.0 * x);
		if (x > x1)
			incidentFlux = 0.0;

		incidentFluxVec.push_back(incidentFlux);
	}

}

