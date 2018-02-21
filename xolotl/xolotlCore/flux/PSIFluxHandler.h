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
	virtual void initializeFluxHandler(IReactionNetwork *network,
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

		// Set the flux index corresponding the the single deterium cluster here
		fluxCluster = network->get(dType, 1);
		// Check that the deuterium cluster is present in the network
		if (!fluxCluster) {
			throw std::string(
					"\nThe single deuterium cluster is not present in the network, "
							"cannot use the flux option!");
		}
		dFluxIndex = fluxCluster->getId() - 1;

		// Set the flux index corresponding the the single tritium cluster here
		fluxCluster = network->get(tType, 1);
		// Check that the tritium cluster is present in the network
		if (!fluxCluster) {
			throw std::string(
					"\nThe single tritium cluster is not present in the network, "
							"cannot use the flux option!");
		}
		tFluxIndex = fluxCluster->getId() - 1;

		return;
	}

	/**
	 * This operation computes the flux due to incoming particles at a given grid point.
	 * \see IFluxHandler.h
	 */
	virtual void computeIncidentFlux(double currentTime,
			double *updatedConcOffset, int xi, int surfacePos) {
		// Recompute the flux vector if a time profile is used
		if (useTimeProfile) {
			fluxAmplitude = getProfileAmplitude(currentTime);
			recomputeFluxHandler(surfacePos);
		}

		if (incidentFluxVec.size() == 0) {
			updatedConcOffset[fluxIndex] += fluxAmplitude;
			updatedConcOffset[dFluxIndex] += fluxAmplitude * 0.55;
			updatedConcOffset[tFluxIndex] += fluxAmplitude * 0.45;
			return;
		}

		// Update the concentration array
		updatedConcOffset[fluxIndex] += incidentFluxVec[xi - surfacePos];
		updatedConcOffset[dFluxIndex] += incidentFluxVec[xi - surfacePos] * 0.55;
		updatedConcOffset[tFluxIndex] += incidentFluxVec[xi - surfacePos] * 0.45;

		return;
	}

};
//end class PSIFluxHandler

}

#endif
