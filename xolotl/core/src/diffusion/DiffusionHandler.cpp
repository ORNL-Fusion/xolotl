#include <xolotl/core/diffusion/DiffusionHandler.h>

namespace xolotl
{
namespace core
{
namespace diffusion
{
void
DiffusionHandler::syncDiffusingClusters(network::IReactionNetwork& network)
{
	using HostUnmanaged =
		Kokkos::View<const IdType*, Kokkos::HostSpace, Kokkos::MemoryUnmanaged>;
	auto clusterIds_h =
		HostUnmanaged(diffusingClusters.data(), diffusingClusters.size());
	diffClusterIds = Kokkos::View<IdType*>(
		"Diffusing Cluster Ids", diffusingClusters.size());
	deep_copy(diffClusterIds, clusterIds_h);

	diffClusters = Kokkos::View<DeviceCluster*>(
		Kokkos::ViewAllocateWithoutInitializing("Diffusing Clusters"),
		diffusingClusters.size());
	auto clusters_h = create_mirror_view(diffClusters);
	for (IdType i = 0; i < diffusingClusters.size(); ++i) {
		clusters_h[i] = network.getClusterCommon(
			diffusingClusters[i], plsm::DeviceMemSpace{});
	}
	deep_copy(diffClusters, clusters_h);
}

void
DiffusionHandler::initialize(
	network::IReactionNetwork& network, std::vector<core::RowColPair>& idPairs)
{
	// Clear the index vector
	diffusingClusters.clear();

	// Consider each cluster
	for (std::size_t i = 0; i < network.getNumClusters(); i++) {
		auto cluster = network.getClusterCommon(i);

		// Get its diffusion factor and migration energy
		double diffFactor = cluster.getDiffusionFactor();
		double migration = cluster.getMigrationEnergy();

		// Don't do anything if the diffusion factor is 0.0
		if (util::equal(diffFactor, 0.0) || migration > migrationThreshold)
			continue;

		// Note that cluster is diffusing.
		diffusingClusters.emplace_back(i);

		// Add a matrix entry for this cluster
		idPairs.push_back(core::RowColPair{i, i});
	}

	this->syncDiffusingClusters(network);
}
} // namespace diffusion
} // namespace core
} // namespace xolotl
