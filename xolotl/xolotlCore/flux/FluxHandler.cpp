#include "FluxHandler.h"
#include <xolotlPerf.h>
#include <iostream>
#include <fstream>
#include <cmath>
#include <limits>
#include <mpi.h>

namespace xolotlCore {

FluxHandler::FluxHandler() :
		stepSize(0.0),
		elementarySurfaceSize(0.0),
		heFluence(0.0),
		heFlux(1.0),
		useTimeProfile(false),
		normFactor(0.0){

}

void FluxHandler::initializeFluxHandler(int nx, double hx, double hy,
		double hz) {

	// Set the step and elementary surface sizes
	stepSize = hx;
	elementarySurfaceSize = hy * hz;

	normFactor = 0.0;
	for (int i = 1; i < nx - 1; i++) {
		double x = (double) i * stepSize;

		normFactor += FitFunction(x) * stepSize;
	}

	// Factor the incident flux will be multiplied by
	double heFluxNormalized = elementarySurfaceSize * heFlux / normFactor;

	// The first value should always be 0.0 because of boundary conditions
	incidentFluxVec.push_back(0.0);

	// Starts a i = 1 because the first value was already put in the vector
	for (int i = 1; i < nx - 1; i++) {
		auto x = i * stepSize;

		auto incidentFlux = heFluxNormalized * FitFunction(x);

		incidentFluxVec.push_back(incidentFlux);
	}

	// The last value should always be 0.0 because of boundary conditions
	incidentFluxVec.push_back(0.0);

	return;
}

void FluxHandler::recomputeFluxHandler() {
	// Factor the incident flux will be multiplied by
	double heFluxNormalized = elementarySurfaceSize * heFlux / normFactor;

	// Get the number of grid points
	int numGridPoints = incidentFluxVec.size();

	// Clear the flux vector
	incidentFluxVec.clear();

	// The first value should always be 0.0 because of boundary conditions
	incidentFluxVec.push_back(0.0);

	// Starts a i = 1 because the first value was already put in the vector
	for (int i = 1; i < numGridPoints; i++) {
		auto x = i * stepSize;

		auto incidentFlux = heFluxNormalized * FitFunction(x);

		incidentFluxVec.push_back(incidentFlux);
	}

	return;
}

void FluxHandler::initializeTimeProfile(std::string fileName) {
	// Set use time profile to true
	useTimeProfile = true;

	// Open file dataFile.dat containing the time and amplitude
	std::ifstream inputFile(fileName.c_str());
	std::string line;

	while (getline(inputFile, line)) {
		if (!line.length() || line[0] == '#')
			continue;
		double xamp = 0.0, yamp = 0.0;
		sscanf(line.c_str(), "%lf %lf", &xamp, &yamp);
		time.push_back(xamp);
		amplitude.push_back(yamp);
	}

	return;
}

double FluxHandler::getAmplitude(double currentTime) const {

	double f = 0.0;

	// if x is smaller than or equal to xi[0]
	if (currentTime <= time[0])
		return f = amplitude[0];

	// if x is greater than or equal to xi[n-1]
	if (currentTime >= time[time.size() - 1])
		return f = amplitude[time.size() - 1];

	// loop to determine the interval x falls in, ie x[k] < x < x[k+1]
	for (int k = 0; k < time.size() - 1; k++) {
		if (currentTime < time[k]) continue;
		if (currentTime > time[k + 1]) continue;

		f = amplitude[k]
				+ (amplitude[k + 1] - amplitude[k]) * (currentTime - time[k])
						/ (time[k + 1] - time[k]);
		break;
	}

	return f;
}

double FluxHandler::getIncidentFlux(std::vector<int> compositionVec,
		std::vector<double> position, double currentTime) {

	// Recompute the flux vector if a time profile is used
	if (useTimeProfile) {
		heFlux = getAmplitude(currentTime);
		recomputeFluxHandler();
	}

	// Get the index number from the position
	int i = position[0] / stepSize;

	// Return the corresponding value
	return incidentFluxVec[i];
}

std::vector<double> FluxHandler::getIncidentFluxVec(double currentTime) {

	// Recompute the flux vector if a time profile is used
	if (useTimeProfile) {
		heFlux = getAmplitude(currentTime);
		recomputeFluxHandler();
	}

	return incidentFluxVec;
}

void FluxHandler::setOutgoingFlux(std::vector<int> compositionVec,
		std::vector<int> position, double time, double outgoingFlux) {

	return;
}

void FluxHandler::incrementHeFluence(double dt) {
	// the fluence is the flux times the time
	heFluence += heFlux * dt;

	return;
}

double FluxHandler::getHeFluence() const {
	return heFluence;
}

void FluxHandler::setHeFlux(double flux) {
	heFlux = flux;
}

double FluxHandler::getHeFlux() const {
	return heFlux;
}

} // end namespace xolotlCore
