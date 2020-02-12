#ifndef FUELFITFLUXHANDLER_H
#define FUELFITFLUXHANDLER_H

#include "FluxHandler.h"
#include <cmath>
#include <MathUtils.h>

namespace xolotlCore {

/**
 * This class realizes the IFluxHandler interface to calculate the incident xenon flux
 * for nuclear fuel.
 */
class FuelFitFluxHandler: public FluxHandler {
private:

	/**
	 * Function that calculate the flux at a given position x (in nm).
	 * This function is not normalized. The surface is supposed to be (100).
	 *
	 * @param x The position where to evaluate he fit
	 * @return The evaluated value
	 */
	double FitFunction(double x) {
		return 1.0;
	}

public:

	/**
	 * The constructor
	 */
	FuelFitFluxHandler() {
	}

	/**
	 * The Destructor
	 */
	~FuelFitFluxHandler() {
	}

	/**
	 * Compute and store the incident flux values at each grid point.
	 * \see IFluxHandler.h
	 */
	void initializeFluxHandler(const IReactionNetwork& network, int surfacePos,
			std::vector<double> grid) {
		// Set the grid
		xGrid = grid;

		// Skip if the flux amplitude is 0.0 and we are not using a time profile
		if (equal(fluxAmplitude, 0.0) && !useTimeProfile) return;

		// Set the flux index corresponding the the single xenon cluster here
		auto fluxCluster = network.get(Species::Xe, 1);
		// Check that the helium cluster is present in the network
		if (!fluxCluster) {
			throw std::string(
					"\nThe single xenon cluster is not present in the network, "
							"cannot use the flux option!");
		}
		fluxIndices.push_back(fluxCluster->getId() - 1);

		return;
	}

	/**
	 * This operation computes the flux due to incoming particles at a given grid point.
	 * \see IFluxHandler.h
	 */
	void computeIncidentFlux(double currentTime,
			double *updatedConcOffset, int xi, int surfacePos) {
		// Skip if no index was set
		if (fluxIndices.size() == 0) return;

		// 0D Case
		if (xGrid.size() == 0) {
			updatedConcOffset[fluxIndices[0]] += fluxAmplitude;
			return;
		}

		// Update the concentration array
		updatedConcOffset[fluxIndices[0]] += fluxAmplitude;

		return;
	}

};
//end class FuelFitFluxHandler

}

#endif
