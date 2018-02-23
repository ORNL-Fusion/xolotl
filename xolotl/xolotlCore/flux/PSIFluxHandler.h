#ifndef PSIFLUXHANDLER_H
#define PSIFLUXHANDLER_H

#include "FluxHandler.h"
#include <vector>
#include <memory>
#include <Constants.h>

namespace xolotlCore {

/**
 * Realizations of this interface are responsible for handling the incident (incoming)
 * flux calculations in tungsten.
 */
class PSIFluxHandler: public FluxHandler {

protected:

	/**
	 * The index of the deuterium cluster.
	 */
	int dFluxIndex;

	/**
	 * The index of the tritium cluster.
	 */
	int tFluxIndex;

public:

	PSIFluxHandler() : dFluxIndex(-1), tFluxIndex(-1) {
	}

	~PSIFluxHandler() {
	}

	/**
	 * Compute and store the incident flux values at each grid point.
	 * \see IFluxHandler.h
	 */
	void initializeFluxHandler(const IReactionNetwork& network, int surfacePos,
			std::vector<double> grid) {
		// Call the general method
		FluxHandler::initializeFluxHandler(network, surfacePos, grid);

		// Set the flux index corresponding the the single helium cluster here
		auto fluxCluster = network.get(Species::He, 1);
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
//end class PSIFluxHandler

}

#endif
