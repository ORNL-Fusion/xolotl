#pragma once

#include <cmath>

#include <xolotl/core/flux/FluxHandler.h>
#include <xolotl/core/network/NEReactionNetwork.h>
#include <xolotl/util/MathUtils.h>

namespace xolotl
{
namespace core
{
namespace flux
{
/**
 * This class realizes the FluxHandler interface to calculate the incident
 * xenon flux for nuclear fuel.
 */
class FuelFitFluxHandler : public FluxHandler
{
private:
	/**
	 * \see FluxHandler.h
	 */
	double
	FitFunction(double x)
	{
		// Constant flux
		return 1.0;
	}

public:
	/**
	 * The constructor
	 */
	FuelFitFluxHandler(const options::IOptions& options) : FluxHandler(options)
	{
		// Change the flux amplitude because we have to take into account
		// that there are one xenon created every 4 fissions.
		this->setFluxAmplitude(
			options.getFluxAmplitude() * options.getFissionYield());
	}

	/**
	 * The Destructor
	 */
	~FuelFitFluxHandler()
	{
	}

	/**
	 * \see IFluxHandler.h
	 */
	void
	initializeFluxHandler(network::IReactionNetwork& network, int surfacePos,
		std::vector<double> grid)
	{
		// Set the grid
		xGrid = grid;

		// Skip if the flux amplitude is 0.0 and we are not using a time profile
		if (util::equal(fluxAmplitude, 0.0) && !useTimeProfile)
			return;

		// Set the flux index corresponding the the single xenon cluster here
		using NetworkType = network::NEReactionNetwork;
		auto& neNetwork = dynamic_cast<NetworkType&>(network);
		NetworkType::Composition comp = NetworkType::Composition::zero();
		comp[NetworkType::Species::Xe] = 1;
		auto cluster = neNetwork.findCluster(comp, plsm::onHost);
		// Check that the helium cluster is present in the network
		if (cluster.getId() == NetworkType::invalidIndex()) {
			throw std::string(
				"\nThe single xenon cluster is not present in the network, "
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
		// Skip if no index was set
		if (fluxIndices.size() == 0)
			return;

		// Update the concentration array
		updatedConcOffset[fluxIndices[0]] += fluxAmplitude;

		return;
	}
};
// end class FuelFitFluxHandler

} // namespace flux
} // namespace core
} // namespace xolotl
