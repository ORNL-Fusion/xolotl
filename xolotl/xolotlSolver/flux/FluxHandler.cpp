#include <FluxHandler.h>
#include "xolotlPerf/xolotlPerf.h"
#include <iostream>

using namespace xolotlSolver;

FluxHandler::FluxHandler() :
		stepSize(0.0e-16), heFluence(0.0e-16), usingMaxHeFluence(false), maxHeFluence(
				0.0e-16), heFlux(1.0), incidentFluxZero(false) {

}

double FluxHandler::getIncidentFlux(std::vector<int> compositionVec,
		std::vector<double> position, double currentTime) {

	auto i = position[0] / stepSize;
	auto incidentFlux = incidentFluxVec[i];

	return incidentFlux;
}

std::vector<double> FluxHandler::getIncidentFluxVec() {
	return incidentFluxVec;
}

void FluxHandler::setOutgoingFlux(std::vector<int> compositionVec,
		std::vector<int> position, double time, double outgoingFlux) {

}

double FluxHandler::incrementHeFluence(double dt, double step) {

	if ( !usingMaxHeFluence || !(heFluence >= maxHeFluence) )
	{
		for (int i = 0; i < incidentFluxVec.size(); i++) {
			heFluence += incidentFluxVec[i] * dt * step;
		}
	}

	if ( !incidentFluxZero && usingMaxHeFluence && heFluence >= maxHeFluence ) {
		for (int i = 0; i < incidentFluxVec.size(); i++)
			incidentFluxVec[i] = 0.0;
		incidentFluxZero = true;
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

bool FluxHandler::getUsingMaxHeFluence() {
	return usingMaxHeFluence;
}

void FluxHandler::setHeFlux(double flux) {
	heFlux = flux;
}

double FluxHandler::getHeFlux() const {
	return heFlux;
}
