#ifndef W111FITFLUXHANDLER_H
#define W111FITFLUXHANDLER_H

#include "FluxHandler.h"
#include <cmath>

namespace xolotlCore {

/**
 * This class realizes the IFluxHandler interface to calculate the incident helium flux
 * for a (111) oriented tungsten material.
 */
class W111FitFluxHandler: public FluxHandler {
private:

	/**
	 * Function that calculate the flux at a given position x (in nm).
	 * This function is not normalized. The surface is supposed to be (111).
	 *
	 * @param x The position where to evaluate he fit
	 * @return The evaluated value
	 */
	double FitFunction(double x) {
		// Value at which the flux goes to 0
		double x1 = 10.0;

		if (x > x1) return 0.0;

		// Compute the fit
		double value = 1.563416 + 7.2071044 * x - 5.4632628 * pow(x, 2)
		+ 1.727342 * pow(x, 3) - 0.3014105 * pow(x, 4) + 0.0311738 * pow(x, 5)
		- 0.0019016 * pow(x, 6) + 0.00006318 * pow(x, 7) - 0.0000008813 * pow(x, 8);

		return value;
	}

public:

	/**
	 * The constructor
	 */
	W111FitFluxHandler() {}

	/**
	 * The Destructor
	 */
	~W111FitFluxHandler() {}

	/**
	 * Compute and store the incident flux values at each grid point.
	 * \see IFluxHandler.h
	 */
	void initializeFluxHandler(IReactionNetwork *network,
			int surfacePos, std::vector<double> grid) {
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
		fluxIndex = fluxCluster->getId() - 1;

		return;
	}

};
//end class W111FitFluxHandler

}

#endif
