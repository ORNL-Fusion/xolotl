#include "FluxHandler.h"
#include <iostream>

using namespace xolotlSolver;

double FluxHandler::getIncidentFlux(std::vector<int> compositionVec,
				std::vector<double> position, double currentTime){

	auto i = position[0] / stepSize;

	return incidentFluxVec[i];
}

void FluxHandler::setOutgoingFlux(std::vector<int> compositionVec,
			std::vector<int> position, double time, double outgoingFlux){

}
