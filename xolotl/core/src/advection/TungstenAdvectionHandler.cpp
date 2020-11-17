#include <xolotl/core/advection/TungstenAdvectionHandler.h>

namespace xolotl
{
namespace core
{
namespace advection
{
TungstenAdvectionHandler::TungstenAdvectionHandler(
	const std::array<double, 7>& sinkStrength) :
	_sinkStrength(sinkStrength)
{
}

TungstenAdvectionHandler::~TungstenAdvectionHandler() = default;

void
TungstenAdvectionHandler::initialize(network::IReactionNetwork& network,
	network::IReactionNetwork::SparseFillMap& ofillMap)
{
	// Clear the index and sink strength vectors
	advectingClusters.clear();
	sinkStrengthVector.clear();

	using NetworkType = network::IPSIReactionNetwork;
	using AmountType = NetworkType::AmountType;

	auto psiNetwork = dynamic_cast<NetworkType*>(&network);
	auto numSpecies = psiNetwork->getSpeciesListSize();
	auto specIdHe = psiNetwork->getHeliumSpeciesId();

	// Initialize the composition
	auto comp = std::vector<AmountType>(numSpecies, 0);

	// Loop on helium clusters from size 1 to 7
	for (std::size_t i = 1; i <= 7; i++) {
		comp[specIdHe()] = i;
		auto clusterId = psiNetwork->findClusterId(comp);
		// Check that the helium cluster is present in the network
		if (clusterId == NetworkType::invalidIndex()) {
			throw std::runtime_error("\nThe helium cluster of size " +
				std::to_string(i) +
				"is not present in the network, "
				"cannot use the advection option!");
		}

		// Get its diffusion coefficient
		double diffFactor =
			psiNetwork->getClusterCommon(clusterId).getDiffusionFactor();

		// Don't do anything if the diffusion factor is 0.0
		if (util::equal(diffFactor, 0.0))
			continue;

		// Switch on the size to get the sink strength (in eV.nm3)
		double sinkStrength = _sinkStrength[i - 1];

		// If the sink strength is still 0.0, this cluster is not advecting
		if (util::equal(sinkStrength, 0.0)) {
			continue;
		}

		// Add it to our collection of advecting clusters.
		advectingClusters.emplace_back(clusterId);

		// Add the sink strength to the vector
		sinkStrengthVector.push_back(sinkStrength);

		// Set the off-diagonal part for the Jacobian to 1
		// Set the ofill value to 1 for this cluster
		ofillMap[clusterId].emplace_back(clusterId);
	}
}
} // namespace advection
} // namespace core
} // namespace xolotl
