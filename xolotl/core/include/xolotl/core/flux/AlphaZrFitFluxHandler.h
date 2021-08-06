#ifndef ALPHAZRFITFLUXHANDLER_H
#define ALPHAZRFITFLUXHANDLER_H

#include <cmath>

#include <xolotl/core/flux/FluxHandler.h>
#include <xolotl/core/network/ZrReactionNetwork.h>

namespace xolotl
{
namespace core
{
namespace flux
{
/**
 * This class realizes the IFluxHandler interface to calculate the incident
 * fluxes for an alpha Zr material.
 */
class AlphaZrFitFluxHandler : public FluxHandler
{
private:
	/**
	 * \see FluxHandler.h
	 */
	double
	FitFunction(double x)
	{
		// Not actually used
		return 0.0;
	}

public:
	/**
	 * The constructor
	 */
	AlphaZrFitFluxHandler(const options::IOptions& options) :
		FluxHandler(options)
	{
	}

	/**
	 * The Destructor
	 */
	~AlphaZrFitFluxHandler()
	{
	}

	/**
	 * \see IFluxHandler.h
	 */
	void
	initializeFluxHandler(network::IReactionNetwork& network, int surfacePos,
		std::vector<double> grid)
	{
		// Only defined in 0D
		if (xGrid.size() == 0) {
			// Add an empty vector
			std::vector<double> tempVector;
			incidentFluxVec.push_back(tempVector);
		}

		// TODO: which clusters have a generation term? Left V_1 and I_1 as
		// examples

		using NetworkType = network::ZrReactionNetwork;
		auto zrNetwork = dynamic_cast<NetworkType*>(&network);

		// Set the flux index corresponding the the single interstitial cluster
		// here
		NetworkType::Composition comp = NetworkType::Composition::zero();
		comp[NetworkType::Species::I] = 1;
		auto cluster = zrNetwork->findCluster(comp, plsm::onHost);
		if (cluster.getId() == NetworkType::invalidIndex()) {
			throw std::runtime_error("\nThe single interstitial cluster is not "
									 "present in the network, "
									 "cannot use the flux option!");
		}
		fluxIndices.push_back(cluster.getId());

		// Look for vacancy now
		comp[NetworkType::Species::I] = 0;
		comp[NetworkType::Species::V] = 1;
		cluster = zrNetwork->findCluster(comp, plsm::onHost);
		if (cluster.getId() == NetworkType::invalidIndex()) {
			throw std::runtime_error(
				"\nThe single vacancy cluster is not present in the network, "
				"cannot use the flux option!");
		}
		fluxIndices.push_back(cluster.getId());

		return;
	}

	/**
	 * \see IFluxHandler.h
	 */
	void
	computeIncidentFlux(
		double currentTime, double* updatedConcOffset, int xi, int surfacePos)
	{
		// Define only for a 0D case
		if (incidentFluxVec[0].size() == 0) {
			updatedConcOffset[fluxIndices[0]] += 0.0; // I1
			updatedConcOffset[fluxIndices[1]] += 0.0; // V1
		}

		else {
			throw std::runtime_error(
				"\nThe alpha Zr problem is not defined for more than 0D!");
		}

		return;
	}
};
// end class AlphaZrFitFluxHandler

} // namespace flux
} // namespace core
} // namespace xolotl

#endif
