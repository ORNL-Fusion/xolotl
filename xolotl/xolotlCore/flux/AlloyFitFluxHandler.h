#ifndef ALLOYFITFLUXHANDLER_H
#define ALLOYFITFLUXHANDLER_H

#include "FluxHandler.h"
#include <cmath>

namespace xolotlCore {

/**
 * This class realizes the IFluxHandler interface to calculate the incident fluxes
 * for the alloy case.
 */
class AlloyFitFluxHandler: public FluxHandler {
private:

	/**
	 * Function that calculate the flux at a given position x (in nm).
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
	AlloyFitFluxHandler() {}

	/**
	 * The Destructor
	 */
	~AlloyFitFluxHandler() {}

	/**
	 * Compute and store the incident flux values at each grid point.
	 * \see IFluxHandler.h
	 */
	void initializeFluxHandler(IReactionNetwork *network,
			int surfacePos, std::vector<double> grid) {
		// Set the grid
		xGrid = grid;

		if (xGrid.size() == 0) return;

		// Clear the flux vector
		incidentFluxVec.clear();
		// The first value corresponding to the surface position should always be 0.0
		incidentFluxVec.push_back(0.0);

		// Starts a i = surfacePos + 1 because the first value was already put in the vector
		for (int i = surfacePos + 1; i < xGrid.size() - 1; i++) {
			// Get the x position
			auto x = xGrid[i] - xGrid[surfacePos];

			// Compute the flux value
			double incidentFlux = FitFunction(x);
			// Add it to the vector
			incidentFluxVec.push_back(incidentFlux);
		}

		// The last value should always be 0.0 because of boundary conditions
		incidentFluxVec.push_back(0.0);

		// Set the flux index corresponding the the single interstitial cluster here
		auto fluxCluster = network->get(iType, 1);
		// Check that the helium cluster is present in the network
		if (!fluxCluster) {
			throw std::string(
					"\nThe single interstitial cluster is not present in the network, "
					"cannot use the flux option!");
		}
		fluxIndex = fluxCluster->getId() - 1;

		return;
	}

	/**
	 * This operation computes the flux due to incoming particles at a given grid point.
	 * \see IFluxHandler.h
	 */
	void computeIncidentFlux(double currentTime, double *updatedConcOffset, int xi, int surfacePos) {
		if (incidentFluxVec.size() == 0) {
			updatedConcOffset[fluxIndex] += fluxAmplitude;
			return;
		}

		// Update the concentration array
		updatedConcOffset[fluxIndex] += incidentFluxVec[xi - surfacePos];

		return;
	}

};
//end class AlloyFitFluxHandler

}

#endif
