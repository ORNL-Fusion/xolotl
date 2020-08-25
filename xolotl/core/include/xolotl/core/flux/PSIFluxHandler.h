#ifndef PSIFLUXHANDLER_H
#define PSIFLUXHANDLER_H

#include <memory>
#include <vector>

#include <xolotl/core/Constants.h>
#include <xolotl/core/flux/FluxHandler.h>
#include <xolotl/core/network/PSIReactionNetwork.h>
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
	 * Compute and store the incident flux values at each grid point.
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

		// Set the flux index corresponding the the single helium cluster here
		using NetworkType =
			network::PSIReactionNetwork<network::PSIFullSpeciesList>;
		auto psiNetwork = dynamic_cast<NetworkType*>(&network);

		// Set the flux index corresponding the the single helium cluster here
		NetworkType::Composition comp;
		// Initialize the composition
		for (auto i : psiNetwork->getSpeciesRange()) {
			comp[i] = 0;
		}
		comp[NetworkType::Species::He] = 1;
		auto cluster = psiNetwork->findCluster(comp, plsm::onHost);
		// Check that the helium cluster is present in the network
		if (cluster.getId() == NetworkType::invalidIndex()) {
			throw std::string(
				"\nThe single helium cluster is not present in the network, "
				"cannot use the flux option!");
		}
		fluxIndices.push_back(cluster.getId());

		return;
	}
};
// end class PSIFluxHandler

} // namespace flux
} // namespace core
} // namespace xolotl

#endif
