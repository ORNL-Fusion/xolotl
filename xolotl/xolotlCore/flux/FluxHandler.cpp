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
	return;
}

void FluxHandler::initializeFluxHandler(int nx, double hx, double hy,
		double hz) {
	// Set the step and elementary surface sizes
	stepSize = hx;
	elementarySurfaceSize = hy * hz;

	// Compute the norm factor because the fit function has an
	// arbitrary amplitude
	normFactor = 0.0;
	// Loop on the x grid points skipping the first and last because
	// of the boundary conditions
	for (int i = 1; i < nx - 1; i++) {
		// Get the x position
		double x = (double) i * stepSize;

		// Add the the value of the function times the step size
		normFactor += FitFunction(x) * stepSize;
	}

	// Factor the incident flux will be multiplied by to get
	// the wanted intensity
	double heFluxNormalized = elementarySurfaceSize * heFlux / normFactor;

	// The first value should always be 0.0 because of boundary conditions
	incidentFluxVec.push_back(0.0);

	// Starts a i = 1 because the first value was already put in the vector
	for (int i = 1; i < nx - 1; i++) {
		// Get the x position
		double x = i * stepSize;

		// Compute the flux value
		double incidentFlux = heFluxNormalized * FitFunction(x);
		// Add it to the vector
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
		// Get the x position
		double x = i * stepSize;

		// Compute the flux value
		double incidentFlux = heFluxNormalized * FitFunction(x);
		// Add it to the vector
		incidentFluxVec.push_back(incidentFlux);
	}

	return;
}

void FluxHandler::initializeTimeProfile(const std::string& fileName) {
	// Set use time profile to true
	useTimeProfile = true;

	// Open file dataFile.dat containing the time and amplitude
	std::ifstream inputFile(fileName.c_str());
	std::string line;

	// Read the file and store the values in the two vectors
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
	// Initialize the amplitude to return
	double f = 0.0;

	// If the time is smaller than or equal than the first stored time
	if (currentTime <= time[0])
		return f = amplitude[0];

	// If the time is larger or equal to the last stored time
	if (currentTime >= time[time.size() - 1])
		return f = amplitude[time.size() - 1];

	// Else loop to determine the interval the time falls in
	// i.e. time[k] < time < time[k + 1]
	for (int k = 0; k < time.size() - 1; k++) {
		if (currentTime < time[k]) continue;
		if (currentTime > time[k + 1]) continue;

		// Compute the amplitude following a linear interpolation between
		// the two stored values
		f = amplitude[k]
				+ (amplitude[k + 1] - amplitude[k]) * (currentTime - time[k])
						/ (time[k + 1] - time[k]);
		break;
	}

	return f;
}

std::vector<double> FluxHandler::getIncidentFluxVec(double currentTime) {
	// Recompute the flux vector if a time profile is used
	if (useTimeProfile) {
		heFlux = getAmplitude(currentTime);
		recomputeFluxHandler();
	}

	return incidentFluxVec;
}

void FluxHandler::incrementHeFluence(double dt) {
	// The fluence is the flux times the time
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
