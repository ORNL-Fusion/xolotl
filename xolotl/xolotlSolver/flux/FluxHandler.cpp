#include <FluxHandler.h>
#include <HandlerRegistryFactory.h>
#include <iostream>

using namespace xolotlSolver;

double FluxHandler::getIncidentFlux(std::vector<int> compositionVec,
		std::vector<double> position, double currentTime) {

	auto i = position[0] / stepSize;
	auto incidentFlux = incidentFluxVec[i];

	if (usingMaxHeFluence && heFluence >= maxHeFluence) {
		incidentFlux = 0.0;
	}

	return incidentFlux;
}

void FluxHandler::setOutgoingFlux(std::vector<int> compositionVec,
		std::vector<int> position, double time, double outgoingFlux) {

}

double FluxHandler::incrementHeFluence(double dt, double step) {

	if (usingMaxHeFluence && heFluence >= maxHeFluence)
		heFluence += 0.0 * dt * step;
	else {
		for (int i = 0; i < incidentFluxVec.size(); i++) {
			heFluence += incidentFluxVec[i] * dt * step;
		}
	}

	return heFluence;
}

double FluxHandler::getHeFluence() const {
	return heFluence;
}

void FluxHandler::setMaxHeFluence(double fluence) {

	usingMaxHeFluence = true;

	maxHeFluence = fluence;
}

double FluxHandler::getMaxHeFluence() const {
	return maxHeFluence;
}

bool FluxHandler::useMaximumHeFluence() {
	return usingMaxHeFluence;
}
