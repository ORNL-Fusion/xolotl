#ifndef PSIFLUXHANDLER_H
#define PSIFLUXHANDLER_H

#include <memory>
#include <vector>

#include <xolotl/core/Constants.h>
#include <xolotl/core/flux/FluxHandler.h>
#include <xolotl/core/network/IPSIReactionNetwork.h>
#include <xolotl/util/MathUtils.h>

namespace xolotl
{
namespace core
{
namespace flux
{
/**
 * Realizations of this interface are responsible for handling the incident
 * (incoming) flux calculations in tungsten.
 */
class PSIFluxHandler : public FluxHandler
{
public:
	PSIFluxHandler(const options::IOptions& options) : FluxHandler(options)
	{
	}

	~PSIFluxHandler()
	{
	}

	/**
	 * \see IFluxHandler.h
	 */
	void
	initializeFluxHandler(network::IReactionNetwork& network, int surfacePos,
		std::vector<double> grid)
	{
		// Call the general method
		FluxHandler::initializeFluxHandler(network, surfacePos, grid);

		// Skip if the flux amplitude is 0.0 and we are not using a time profile
		if (util::equal(fluxAmplitude, 0.0) && !useTimeProfile)
			return;

		// Set the flux index corresponding the single helium cluster here
		using NetworkType = network::IPSIReactionNetwork;
		using AmountType = NetworkType::AmountType;

		auto psiNetwork = dynamic_cast<NetworkType*>(&network);
		auto numSpecies = psiNetwork->getSpeciesListSize();
		auto specIdHe = psiNetwork->getHeliumSpeciesId();

		auto comp = std::vector<AmountType>(numSpecies, 0);
		comp[specIdHe()] = 1;
		auto clusterId = psiNetwork->findClusterId(comp);
		// Check that the helium cluster is present in the network
		if (clusterId == NetworkType::invalidIndex()) {
			throw std::runtime_error(
				"\nThe single helium cluster is not present in the network, "
				"cannot use the flux option!");
		}
		fluxIndices.push_back(clusterId);
	}
};
// end class PSIFluxHandler

} // namespace flux
} // namespace core
} // namespace xolotl

#endif
