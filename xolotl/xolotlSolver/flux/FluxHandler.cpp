#include <FluxHandler.h>
#include <HandlerRegistryFactory.h>
#include <iostream>
#include <limits>

using namespace xolotlSolver;

FluxHandler::FluxHandler() :
		stepSize(0.0),
		heFluence(0.0),
		usingMaxHeFluence(false),
		maxHeFluence(std::numeric_limits<double>::max()),
		heFlux(1.0),
		incidentFluxZero(false) {

}

double FluxHandler::getIncidentFlux(std::vector<int> compositionVec,
		std::vector<double> position, double currentTime) {

	// Get the index number from the position
	int i = position[0] / stepSize;

	// Return the corresponding value
	return incidentFluxVec[i];
}

std::vector<double> FluxHandler::getIncidentFluxVec() {
	return incidentFluxVec;
}

void FluxHandler::setOutgoingFlux(std::vector<int> compositionVec,
		std::vector<int> position, double time, double outgoingFlux) {

	return;
}

void FluxHandler::incrementHeFluence(double dt) {

	if (heFluence < maxHeFluence)
	{
		heFluence += heFlux * dt;
	}

	else if (!incidentFluxZero) {
		for (int i = 0; i < incidentFluxVec.size(); i++)
			incidentFluxVec[i] = 0.0;
		incidentFluxZero = true;
	}

	return;
}

double FluxHandler::getHeFluence() const {
	return heFluence;
}

void FluxHandler::setMaxHeFluence(double fluence) {
	usingMaxHeFluence = true;
	maxHeFluence = fluence;

	return;
}

double FluxHandler::getMaxHeFluence() const {
	return maxHeFluence;
}

bool FluxHandler::getUsingMaxHeFluence() {
	return usingMaxHeFluence;
}

void FluxHandler::setHeFlux(double flux) {
	heFlux = flux;
}

double FluxHandler::getHeFlux() const {
	return heFlux;
}
