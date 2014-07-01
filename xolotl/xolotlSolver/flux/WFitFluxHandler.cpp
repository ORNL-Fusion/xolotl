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
		auto incidentFlux = 0.0006 * x * x * x - 0.0087 * x * x + 0.0300 * x;
		incidentFluxVec.push_back(incidentFlux);
	}
//	std::cout << "\n\nincidentFluxVec: " << std::endl;
//	for(int i = 0; i < numGridpoints; i++)
//		std::cout << incidentFluxVec[i] << std::endl;

}

