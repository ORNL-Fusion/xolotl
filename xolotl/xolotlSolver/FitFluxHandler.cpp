#include "FitFluxHandler.h"

using namespace xolotlSolver;

FitFluxHandler::FitFluxHandler(){

}

FitFluxHandler::~FitFluxHandler() {

}

double FitFluxHandler::getIncidentFlux(std::vector<int> compositionVec,
				std::vector<double> position, double currentTime){

	double incidentFlux = 0.0e-16;

	incidentFlux = 0.0006 * position[1] * position[1] * position[1] - 0.0087 * position[1] * position[1] + 0.0300 * position[1];

	return incidentFlux;
}

void FitFluxHandler::setOutgoingFlux(std::vector<int> compositionVec,
			std::vector<int> position, double time, double outgoingFlux){

}
