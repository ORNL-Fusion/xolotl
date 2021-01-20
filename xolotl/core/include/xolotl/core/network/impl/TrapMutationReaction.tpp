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
	Superclass(reactionData, clusterData, reactionId),
	_heClId(cluster0),
	_heVClId(cluster1),
	_iClId(cluster2)
{
	this->initialize();
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
	return 0.0;
}

template <typename TNetwork, typename TDerived>
KOKKOS_INLINE_FUNCTION
double
TrapMutationReaction<TNetwork, TDerived>::computeRate(double largestRate)
{
	const auto& desorp = this->_clusterData.desorption();
	if (_heClId == desorp.id) {
		return (1.0 - desorp.portion) / desorp.portion;
	}

	// Multiply the biggest rate in the network by 1000.0
	// so that trap-mutation overcomes any other reaction
	return (1000.0 * largestRate);
}

template <typename TNetwork, typename TDerived>
KOKKOS_INLINE_FUNCTION
void
TrapMutationReaction<TNetwork, TDerived>::computeConnectivity(
	const Connectivity& connectivity)
{
	this->addConnectivity(_heClId, _heClId, connectivity);
	this->addConnectivity(_heVClId, _heClId, connectivity);
	this->addConnectivity(_iClId, _heClId, connectivity);
}

template <typename TNetwork, typename TDerived>
KOKKOS_INLINE_FUNCTION
void
TrapMutationReaction<TNetwork, TDerived>::computeReducedConnectivity(
	const Connectivity& connectivity)
{
	this->addConnectivity(_heClId, _heClId, connectivity);
}

template <typename TNetwork, typename TDerived>
KOKKOS_INLINE_FUNCTION
double
TrapMutationReaction<TNetwork, TDerived>::getAppliedRate(
	IndexType gridIndex) const
{
	double rate = this->_rate[gridIndex];
	if (_heClId == this->_clusterData.desorption().id) {
		rate *= this->_clusterData.currentDesorpLeftSideRate();
	}
	return rate;
}

template <typename TNetwork, typename TDerived>
KOKKOS_INLINE_FUNCTION
void
TrapMutationReaction<TNetwork, TDerived>::computeFlux(
	ConcentrationsView concentrations, FluxesView fluxes, IndexType gridIndex)
{
	auto rate = getAppliedRate(gridIndex);
	auto f = rate * concentrations[_heClId] *
		this->_clusterData.currentDisappearingRate();

	Kokkos::atomic_sub(&fluxes[_heClId], f);
	Kokkos::atomic_add(&fluxes[_heVClId], f);
	Kokkos::atomic_add(&fluxes[_iClId], f);
}

template <typename TNetwork, typename TDerived>
KOKKOS_INLINE_FUNCTION
void
TrapMutationReaction<TNetwork, TDerived>::computePartialDerivatives(
	ConcentrationsView concentrations, Kokkos::View<double*> values,
	Connectivity connectivity, IndexType gridIndex)
{
	auto rate = getAppliedRate(gridIndex);

	Kokkos::atomic_sub(&values(connectivity(_heClId, _heClId)), rate);
	Kokkos::atomic_add(&values(connectivity(_heVClId, _heClId)), rate);
	Kokkos::atomic_add(&values(connectivity(_iClId, _heClId)), rate);
}

template <typename TNetwork, typename TDerived>
KOKKOS_INLINE_FUNCTION
void
TrapMutationReaction<TNetwork, TDerived>::computeReducedPartialDerivatives(
	ConcentrationsView concentrations, Kokkos::View<double*> values,
	Connectivity connectivity, IndexType gridIndex)
{
	auto rate = getAppliedRate(gridIndex);
	Kokkos::atomic_sub(&values(connectivity(_heClId, _heClId)), rate);
}

template <typename TNetwork, typename TDerived>
KOKKOS_INLINE_FUNCTION
double
TrapMutationReaction<TNetwork, TDerived>::computeLeftSideRate(
	ConcentrationsView concentrations, IndexType clusterId, IndexType gridIndex)
{
	return 0.0;
}
} // namespace network
} // namespace core
} // namespace xolotl
