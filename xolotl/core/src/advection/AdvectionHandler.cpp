#include <xolotl/core/advection/AdvectionHandler.h>

namespace xolotl
{
namespace core
{
namespace advection
{
void
AdvectionHandler::syncAdvectingClusters(network::IReactionNetwork& network)
{
	using HostUnmanaged =
		Kokkos::View<const IdType*, Kokkos::HostSpace, Kokkos::MemoryUnmanaged>;
	auto clusterIds_h =
		HostUnmanaged(advectingClusters.data(), advectingClusters.size());
	advClusterIds = Kokkos::View<IdType*>(
		Kokkos::ViewAllocateWithoutInitializing("Advecting Cluster Ids"),
		advectingClusters.size());
	deep_copy(advClusterIds, clusterIds_h);

	using DeviceCluster = network::ClusterCommon<plsm::DeviceMemSpace>;
	advClusters = Kokkos::View<DeviceCluster*>(
		Kokkos::ViewAllocateWithoutInitializing("Advecting Clusters"),
		advectingClusters.size());
	auto clusters_h = create_mirror_view(advClusters);
	for (IdType i = 0; i < advectingClusters.size(); ++i) {
		clusters_h[i] = network.getClusterCommon(
			advectingClusters[i], plsm::DeviceMemSpace{});
	}
	deep_copy(advClusters, clusters_h);
}

void
AdvectionHandler::syncSinkStrengths()
{
	auto sinkStrengths_h =
		Kokkos::View<const double*, Kokkos::HostSpace, Kokkos::MemoryUnmanaged>(
			sinkStrengthVector.data(), sinkStrengthVector.size());
	advSinkStrengths =
		Kokkos::View<double*>("Sink Strengths", sinkStrengthVector.size());
	deep_copy(advSinkStrengths, sinkStrengths_h);
}
} // namespace advection
} // namespace core
} // namespace xolotl
