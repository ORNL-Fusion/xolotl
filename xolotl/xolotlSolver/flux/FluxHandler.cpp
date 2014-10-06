#include <FluxHandler.h>
#include <xolotlPerf.h>
#include <iostream>
#include <limits>
#include <mpi.h>

using namespace xolotlSolver;

FluxHandler::FluxHandler() :
		stepSize(0.0),
		heFluence(0.0),
		usingMaxHeFluence(false),
		maxHeFluence(std::numeric_limits<double>::max()),
		heFlux(1.0),
		incidentFluxZero(false) {

}

void FluxHandler::initializeFluxHandler(int numGridpoints, double step) {

	// Set the step size
	stepSize = step;

	double normFactor = 0.0;
	for (int i = 1; i < numGridpoints; i++) {
		double x = (double) i * stepSize;

		normFactor += FitFunction(x) * stepSize;
	}

	// Factor the incident flux will be multiplied by
	double heFluxNormalized = heFlux / normFactor;

	// The first value should always be 0.0 because of boundary conditions
	incidentFluxVec.push_back(0.0);

	// Starts a i = 1 because the first value was already put in the vector
	for (int i = 1; i < numGridpoints; i++) {
		auto x = i * stepSize;

		auto incidentFlux = heFluxNormalized * FitFunction(x);

		incidentFluxVec.push_back(incidentFlux);
	}

	return;
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
