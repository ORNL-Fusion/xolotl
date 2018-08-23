#ifndef FUELFITFLUXHANDLER_H
#define FUELFITFLUXHANDLER_H

#include "FluxHandler.h"
#include <cmath>

namespace xolotlCore {

/**
 * This class realizes the IFluxHandler interface to calculate the incident xenon flux
 * for nuclear fuel.
 */
class FuelFitFluxHandler: public FluxHandler {
private:

	/**
	 * Function that calculate the flux at a given position x (in nm).
	 * This function is not normalized.
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
	void initializeFluxHandler(IReactionNetwork *network, int surfacePos,
			std::vector<double> grid) {
		// Call the general method
		FluxHandler::initializeFluxHandler(network, surfacePos, grid);

		// Set the flux index corresponding the the single xenon cluster here
		auto fluxCluster = network->get(xeType, 1);
		// Check that the helium cluster is present in the network
		if (!fluxCluster) {
			throw std::string(
					"\nThe single xenon cluster is not present in the network, "
							"cannot use the flux option!");
		}
		fluxIndex = fluxCluster->getId() - 1;

		return;
	}

};
//end class FuelFitFluxHandler

}

#endif
