#pragma once

namespace xolotl
{
namespace core
{
namespace network
{
template <typename TNetwork, typename TDerived>
KOKKOS_INLINE_FUNCTION
TrapMutationReaction<TNetwork, TDerived>::TrapMutationReaction(
	ReactionDataRef reactionData, ClusterDataRef clusterData,
	IndexType reactionId, IndexType cluster0, IndexType cluster1,
	IndexType cluster2) :
	Superclass(reactionData, clusterData, reactionId)
{
}

template <typename TNetwork, typename TDerived>
KOKKOS_INLINE_FUNCTION
TrapMutationReaction<TNetwork, TDerived>::TrapMutationReaction(
	ReactionDataRef reactionData, ClusterDataRef clusterData,
	IndexType reactionId, const detail::ClusterSet& clusterSet) :
	TrapMutationReaction(reactionData, clusterData, reactionId,
		clusterSet.cluster0, clusterSet.cluster1, clusterSet.cluster2)
{
}

template <typename TNetwork, typename TDerived>
KOKKOS_INLINE_FUNCTION
double
TrapMutationReaction<TNetwork, TDerived>::computeRate(IndexType gridIndex)
{
	// TODO
	return 0.0;
}

template <typename TNetwork, typename TDerived>
KOKKOS_INLINE_FUNCTION
void
TrapMutationReaction<TNetwork, TDerived>::computeConnectivity(
	const Connectivity& connectivity)
{
	// TODO
}

template <typename TNetwork, typename TDerived>
KOKKOS_INLINE_FUNCTION
void
TrapMutationReaction<TNetwork, TDerived>::computeReducedConnectivity(
	const Connectivity& connectivity)
{
	// TODO
}

template <typename TNetwork, typename TDerived>
KOKKOS_INLINE_FUNCTION
void
TrapMutationReaction<TNetwork, TDerived>::computeFlux(
	ConcentrationsView concentrations, FluxesView fluxes, IndexType gridIndex)
{
	// TODO
}

template <typename TNetwork, typename TDerived>
KOKKOS_INLINE_FUNCTION
void
TrapMutationReaction<TNetwork, TDerived>::computePartialDerivatives(
	ConcentrationsView concentrations, Kokkos::View<double*> values,
	Connectivity connectivity, IndexType gridIndex)
{
	// TODO
}

template <typename TNetwork, typename TDerived>
KOKKOS_INLINE_FUNCTION
void
TrapMutationReaction<TNetwork, TDerived>::computeReducedPartialDerivatives(
	ConcentrationsView concentrations, Kokkos::View<double*> values,
	Connectivity connectivity, IndexType gridIndex)
{
	// TODO
}

template <typename TNetwork, typename TDerived>
KOKKOS_INLINE_FUNCTION
double
TrapMutationReaction<TNetwork, TDerived>::computeLeftSideRate(
	ConcentrationsView concentrations, IndexType clusterId, IndexType gridIndex)
{
	// TODO
    return 0.0;
}
} // namespace network
} // namespace core
} // namespace xolotl
