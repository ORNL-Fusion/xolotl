#include "WFitFluxHandler.h"
#include <iostream>

using namespace xolotlSolver;

WFitFluxHandler::WFitFluxHandler(){

}

void WFitFluxHandler::initializeFluxHandler(int numGridpoints, double step){

	// Set the step size
	stepSize = step;

	for(int i = 0; i < numGridpoints; i++)
	{
		auto x = i * stepSize;
		auto incidentFlux = 6.0 * x * x * x - 87.0 * x * x + 300.0 * x;
		if (x > 7.0) incidentFlux = 0.0;
		
		// Insure the flux is nonnegative
		incidentFlux = std::max(0.0, incidentFlux);

		incidentFluxVec.push_back(incidentFlux);
	}

}

