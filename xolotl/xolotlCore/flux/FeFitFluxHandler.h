#ifndef FEFITFLUXHANDLER_H
#define FEFITFLUXHANDLER_H

#include "FluxHandler.h"
#include <cmath>

namespace xolotlCore {

/**
 * This class realizes the IFluxHandler interface to calculate the incident fluxes
 * for an iron material.
 */
class FeFitFluxHandler: public FluxHandler {
private:

	/**
	 * Function that calculate the flux at a given position x (in nm).
	 * This function is not normalized. This is for iron.
	 *
	 * @param x The position where to evaluate he fit
	 * @return The evaluated value
	 */
	double FitFunction(double x) {
		return 0.0;
	}

public:

	/**
	 * The constructor
	 */
	FeFitFluxHandler() {}

	/**
	 * The Destructor
	 */
	~FeFitFluxHandler() {}

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
		fluxIndices.push_back(fluxCluster->getId() - 1);

		// Look for interstitial now
		fluxCluster = network->get(iType, 1);
		if (!fluxCluster) {
			throw std::string(
					"\nThe single interstitial cluster is not present in the network, "
					"cannot use the flux option!");
		}
		fluxIndices.push_back(fluxCluster->getId() - 1);

		// Look for vacancies now
		fluxCluster = network->get(vType, 1);
		if (!fluxCluster) {
			throw std::string(
					"\nThe single vacancy cluster is not present in the network, "
					"cannot use the flux option!");
		}
		fluxIndices.push_back(fluxCluster->getId() - 1);
		fluxCluster = network->get(vType, 2);
		if (!fluxCluster) {
			throw std::string(
					"\nThe double vacancy cluster is not present in the network, "
					"cannot use the flux option!");
		}
		fluxIndices.push_back(fluxCluster->getId() - 1);
		fluxCluster = network->get(vType, 3);
		if (!fluxCluster) {
			throw std::string(
					"\nThe triple vacancy cluster is not present in the network, "
					"cannot use the flux option!");
		}
		fluxIndices.push_back(fluxCluster->getId() - 1);
		fluxCluster = network->get(vType, 4);
		if (!fluxCluster) {
			throw std::string(
					"\nThe quadruple vacancy cluster is not present in the network, "
					"cannot use the flux option!");
		}
		fluxIndices.push_back(fluxCluster->getId() - 1);
		fluxCluster = network->get(vType, 5);
		if (!fluxCluster) {
			throw std::string(
					"\nVacancy 5 cluster is not present in the network, "
					"cannot use the flux option!");
		}
		fluxIndices.push_back(fluxCluster->getId() - 1);
		fluxCluster = network->get(vType, 9);
		if (!fluxCluster) {
			throw std::string(
					"\nVacancy 9 cluster is not present in the network, "
					"cannot use the flux option!");
		}
		fluxIndices.push_back(fluxCluster->getId() - 1);

		return;
	}

	/**
	 * This operation computes the flux due to incoming particles at a given grid point.
	 * \see IFluxHandler.h
	 */
	void computeIncidentFlux(double currentTime, double *updatedConcOffset, int xi, int surfacePos) {
		// Define only for a 0D case
		if (incidentFluxVec.size() == 0) {
			updatedConcOffset[fluxIndices[0]] += 2.11e-11; // He1
			updatedConcOffset[fluxIndices[1]] += 1.49e-05; // I1
			updatedConcOffset[fluxIndices[2]] += 9.91e-06; // V1
			updatedConcOffset[fluxIndices[3]] += 1.51e-06; // V2
			updatedConcOffset[fluxIndices[4]] += 2.60e-07; // V3
			updatedConcOffset[fluxIndices[5]] += 1.58e-07; // V4
			updatedConcOffset[fluxIndices[6]] += 6.29e-08; // V5
			updatedConcOffset[fluxIndices[7]] += 3.16e-08; // V9
		}

		else {
			throw std::string(
					"\nThe iron problem is not defined for more than 0D!");
		}

		return;
	}

};
//end class FeFitFluxHandler

}

#endif
