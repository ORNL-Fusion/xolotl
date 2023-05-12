#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>

#include <mpi.h>

#include <xolotl/core/flux/FluxHandler.h>

namespace xolotl
{
namespace core
{
namespace flux
{
FluxHandler::FluxHandler(const options::IOptions& options) :
	fluence(0.0),
	fluxAmplitude(0.0),
	useTimeProfile(false),
	normFactor(0.0)
{
	// Wrong if both the flux and time profile options are used
	if (options.useFluxAmplitude() && options.useFluxTimeProfile()) {
		throw std::runtime_error("\nA constant flux value AND a time profile "
								 "cannot both be given.");
	}
	else if (options.useFluxAmplitude()) {
		// Set the constant value of the flux
		setFluxAmplitude(options.getFluxAmplitude());
	}
	else if (options.useFluxTimeProfile()) {
		// Initialize the time profile
		initializeTimeProfile(options.getFluxTimeProfileFilePath());
	}
}

void
FluxHandler::initializeFluxHandler(network::IReactionNetwork& network,
	int surfacePos, std::vector<double> grid)
{
	// Set the grid
	xGrid = grid;

	if (xGrid.size() == 0) {
		// Add an empty vector
		std::vector<double> tempVector;
		incidentFluxVec.push_back(tempVector);

		return;
	}

	// Compute the norm factor because the fit function has an
	// arbitrary amplitude
	normFactor = 0.0;
	// Loop on the x grid points skipping the first after the surface position
	// and last because of the boundary conditions
	for (int i = surfacePos + 1; i < xGrid.size() - 3; i++) {
		// Get the x position
		double x = (xGrid[i] + xGrid[i + 1]) / 2.0 - xGrid[surfacePos + 1];

		// Add the the value of the function times the step size
		normFactor += FitFunction(x) * (xGrid[i + 1] - xGrid[i]);
	}

	// Factor the incident flux will be multiplied by to get
	// the wanted intensity
	double fluxNormalized = 0.0;
	if (normFactor > 0.0)
		fluxNormalized = fluxAmplitude / normFactor;

	// Clear the flux vector
	incidentFluxVec.clear();
	// The first value corresponding to the surface position should always be
	// 0.0
	std::vector<double> tempVector;
	tempVector.push_back(0.0);

	// Starts a i = surfacePos + 1 because the first value was already put in
	// the vector
	for (int i = surfacePos + 1; i < xGrid.size() - 3; i++) {
		// Get the x position
		auto x = (xGrid[i] + xGrid[i + 1]) / 2.0 - xGrid[surfacePos + 1];

		// Compute the flux value
		double incidentFlux = fluxNormalized * FitFunction(x);
		// Add it to the vector
		tempVector.push_back(incidentFlux);
	}

	// The last value should always be 0.0 because of boundary conditions
	tempVector.push_back(0.0);

	// Add it to the vector of fluxes
	incidentFluxVec.push_back(tempVector);
}

void
FluxHandler::recomputeFluxHandler(int surfacePos)
{
	// Factor the incident flux will be multiplied by
	double fluxNormalized = 0.0;
	if (normFactor > 0.0)
		fluxNormalized = fluxAmplitude / normFactor;

	// Starts at i = surfacePos + 1 because the first values were already put in
	// the vector
	for (int i = surfacePos + 1; i < xGrid.size() - 3; i++) {
		// Get the x position
		auto x = (xGrid[i] + xGrid[i + 1]) / 2.0 - xGrid[surfacePos + 1];

		// Compute the flux value
		double incidentFlux = fluxNormalized * FitFunction(x);
		// Add it to the vector
		incidentFluxVec[0][i - surfacePos] = incidentFlux;
	}
}

void
FluxHandler::syncFluxIndices()
{
	auto ids_h =
		Kokkos::View<IdType*, Kokkos::HostSpace, Kokkos::MemoryUnmanaged>(
			fluxIndices.data(), fluxIndices.size());
	fluxIds = Kokkos::View<IdType*>(
		Kokkos::ViewAllocateWithoutInitializing("Flux Indices"),
		fluxIndices.size());
	deep_copy(fluxIds, ids_h);
}

void
FluxHandler::syncIncidentFluxVec()
{
	incidentFlux = Kokkos::View<double**>(
		"Incident Flux Vec", incidentFluxVec.size(), incidentFluxVec[0].size());
	auto incidentFlux_h = create_mirror_view(incidentFlux);
	for (std::size_t i = 0; i < incidentFluxVec.size(); ++i) {
		for (std::size_t j = 0; j < incidentFluxVec[i].size(); ++j) {
			incidentFlux_h(i, j) = incidentFluxVec[i][j];
		}
	}
	deep_copy(incidentFlux, incidentFlux_h);
}

void
FluxHandler::initializeTimeProfile(const std::string& fileName)
{
	// Set use time profile to true
	useTimeProfile = true;

	// Open file containing the time and amplitude
	std::ifstream inputFile(fileName.c_str());
	std::string line;

	// Read the file and store the values in the two vectors
	while (getline(inputFile, line)) {
		if (!line.length() || line[0] == '#')
			continue;
		double xamp = 0.0, yamp = 0.0;
		sscanf(line.c_str(), "%lf %lf", &xamp, &yamp);
		time.push_back(xamp);
		amplitudes.push_back(yamp);
	}
}

double
FluxHandler::getProfileAmplitude(double currentTime) const
{
	// Initialize the amplitude to return
	double f = 0.0;

	// If the time is smaller than or equal than the first stored time
	if (currentTime <= time[0])
		return f = amplitudes[0];

	// If the time is larger or equal to the last stored time
	if (currentTime >= time[time.size() - 1])
		return f = amplitudes[time.size() - 1];

	// Else loop to determine the interval the time falls in
	// i.e. time[k] < time < time[k + 1]
	for (unsigned int k = 0; k < time.size() - 1; k++) {
		if (currentTime < time[k])
			continue;
		if (currentTime > time[k + 1])
			continue;

		// Compute the amplitude following a linear interpolation between
		// the two stored values
		f = amplitudes[k] +
			(amplitudes[k + 1] - amplitudes[k]) * (currentTime - time[k]) /
				(time[k + 1] - time[k]);
		break;
	}

	return f;
}

void
FluxHandler::computeIncidentFlux(double currentTime,
	Kokkos::View<double*> updatedConcOffset, int xi, int surfacePos)
{
	// Skip if no index was set
	if (fluxIndices.size() == 0)
		return;

	// Recompute the flux vector if a time profile is used
	if (useTimeProfile) {
		fluxAmplitude = getProfileAmplitude(currentTime);
		recomputeFluxHandler(surfacePos);
	}

	double value{};
	if (incidentFluxVec[0].size() == 0) {
		value = fluxAmplitude;
	}
	else {
		value = incidentFluxVec[0][xi - surfacePos];
	}

	// Update the concentration array
	auto id = fluxIndices[0];
	Kokkos::parallel_for(
		1, KOKKOS_LAMBDA(std::size_t) { updatedConcOffset[id] += value; });
}

void
FluxHandler::incrementFluence(double dt)
{
	// The fluence is the flux times the time
	fluence += fluxAmplitude * dt;
}

void
FluxHandler::computeFluence(double time)
{
	// The fluence is the flux times the time
	fluence = fluxAmplitude * time;
}

double
FluxHandler::getFluence() const
{
	return fluence;
}

void
FluxHandler::setFluxAmplitude(double flux)
{
	fluxAmplitude = flux;
}

double
FluxHandler::getFluxAmplitude() const
{
	return fluxAmplitude;
}

double
FluxHandler::getFluxRate() const
{
	if (incidentFluxVec[0].size() == 0)
		return fluxAmplitude;
	return fluxAmplitude / normFactor;
}

std::vector<double>
FluxHandler::getInstantFlux(double time) const
{
	// Create the vector to return
	std::vector<double> toReturn;

	// Get the current amplitude
	auto ampl = getProfileAmplitude(time);
	toReturn.push_back(ampl);

	return toReturn;
}

std::vector<IdType>
FluxHandler::getFluxIndices() const
{
	return fluxIndices;
}

} // end namespace flux
} // end namespace core
} // end namespace xolotl
