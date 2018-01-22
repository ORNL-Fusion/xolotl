#ifndef W100FITFLUXHANDLER_H
#define W100FITFLUXHANDLER_H

#include "FluxHandler.h"
#include <cmath>

namespace xolotlCore {

/**
 * This class realizes the IFluxHandler interface to calculate the incident helium flux
 * for a (100) oriented tungsten material.
 */
class W100FitFluxHandler: public FluxHandler {
private:

	/**
	 * Function that calculate the flux at a given position x (in nm).
	 * This function is not normalized. The surface is supposed to be (100).
	 *
	 * @param x The position where to evaluate he fit
	 * @return The evaluated value
	 */
	double FitFunction(double x) {
		// Value at which the flux goes to 0
		double x1 = 10.0;

		if (x > x1)
			return 0.0;

		// Compute the fit
		double value = 7.00876507 + 0.6052078 * x - 3.01711048 * pow(x, 2)
				+ 1.36595786 * pow(x, 3) - 0.295595 * pow(x, 4)
				+ 0.03597462 * pow(x, 5) - 0.0025142 * pow(x, 6)
				+ 0.0000942235 * pow(x, 7) - 0.0000014679 * pow(x, 8);

		return value;
	}

public:

	/**
	 * The constructor
	 */
	W100FitFluxHandler() {
	}

	/**
	 * The Destructor
	 */
	~W100FitFluxHandler() {
	}

	/**
	 * Compute and store the incident flux values at each grid point.
	 * \see IFluxHandler.h
	 */
	void initializeFluxHandler(IReactionNetwork *network, int surfacePos,
			std::vector<double> grid) {
		// Call the general method
		FluxHandler::initializeFluxHandler(network, surfacePos, grid);

		// Set the flux index corresponding the the single helium cluster here
		auto fluxCluster = network->get(heType, 1);
		// Check that the helium cluster is present in the network
		if (!fluxCluster) {
			throw std::string(
					"\nThe single helium cluster is not present in the network, "
							"cannot use the flux option!");
		}
		fluxIndices.push_back(fluxCluster->getId() - 1);

		return;
	}

};
//end class W100FitFluxHandler

}

#endif
