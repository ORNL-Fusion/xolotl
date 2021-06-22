#pragma once

namespace xolotl
{
namespace core
{
namespace network
{
template <typename TNetwork, typename TDerived>
KOKKOS_INLINE_FUNCTION
BurstingReaction<TNetwork, TDerived>::BurstingReaction(
	ReactionDataRef reactionData, ClusterDataRef clusterData,
	IndexType reactionId, IndexType cluster0, IndexType cluster1) :
	Superclass(reactionData, clusterData, reactionId),
	_reactant(cluster0),
	_product(cluster1)
{
	this->initialize();
}

template <typename TNetwork, typename TDerived>
KOKKOS_INLINE_FUNCTION
BurstingReaction<TNetwork, TDerived>::BurstingReaction(
	ReactionDataRef reactionData, ClusterDataRef clusterData,
	IndexType reactionId, const detail::ClusterSet& clusterSet) :
	BurstingReaction(reactionData, clusterData, reactionId, clusterSet.cluster0,
		clusterSet.cluster1)
{
}

template <typename TNetwork, typename TDerived>
KOKKOS_INLINE_FUNCTION
double
BurstingReaction<TNetwork, TDerived>::computeRate(IndexType gridIndex)
{
	return 1.0e9;
}

template <typename TNetwork, typename TDerived>
KOKKOS_INLINE_FUNCTION
void
BurstingReaction<TNetwork, TDerived>::updateRates(double largestRate)
{
	for (IndexType i = 0; i < this->_rate.extent(0); ++i) {
		this->_rate(i) = 1.0e9;
	}
}

template <typename TNetwork, typename TDerived>
KOKKOS_INLINE_FUNCTION
void
BurstingReaction<TNetwork, TDerived>::computeConnectivity(
	const Connectivity& connectivity)
{
	this->addConnectivity(_reactant, _reactant, connectivity);
	if (_product != invalidIndex)
		this->addConnectivity(_product, _reactant, connectivity);
}

template <typename TNetwork, typename TDerived>
KOKKOS_INLINE_FUNCTION
void
BurstingReaction<TNetwork, TDerived>::computeReducedConnectivity(
	const Connectivity& connectivity)
{
	this->addConnectivity(_reactant, _reactant, connectivity);
}

template <typename TNetwork, typename TDerived>
KOKKOS_INLINE_FUNCTION
double
BurstingReaction<TNetwork, TDerived>::getAppliedRate(IndexType gridIndex) const
{
	// Get the radius of the cluster
	auto cl = this->_clusterData.getCluster(_reactant);
	auto radius = cl.getReactionRadius();

	// Get the current depth
	auto depth = this->_clusterData.getDepth();
	auto tau = this->_clusterData.getTauBursting();
	auto f = this->_clusterData.getFBursting();
	return f * (radius / depth) *
		util::min(1.0, exp(-(depth - tau) / (2.0 * tau)));
}

template <typename TNetwork, typename TDerived>
KOKKOS_INLINE_FUNCTION
void
BurstingReaction<TNetwork, TDerived>::computeFlux(
	ConcentrationsView concentrations, FluxesView fluxes, IndexType gridIndex)
{
	auto rate = getAppliedRate(gridIndex);
	auto f = rate * concentrations[_reactant];

	Kokkos::atomic_sub(&fluxes[_reactant], f);
	if (_product != invalidIndex)
		Kokkos::atomic_add(&fluxes[_product], f);
}

template <typename TNetwork, typename TDerived>
KOKKOS_INLINE_FUNCTION
void
BurstingReaction<TNetwork, TDerived>::computePartialDerivatives(
	ConcentrationsView concentrations, Kokkos::View<double*> values,
	Connectivity connectivity, IndexType gridIndex)
{
	auto rate = getAppliedRate(gridIndex);

	Kokkos::atomic_sub(&values(connectivity(_reactant, _reactant)), rate);
	if (_product != invalidIndex)
		Kokkos::atomic_add(&values(connectivity(_product, _reactant)), rate);
}

template <typename TNetwork, typename TDerived>
KOKKOS_INLINE_FUNCTION
void
BurstingReaction<TNetwork, TDerived>::computeReducedPartialDerivatives(
	ConcentrationsView concentrations, Kokkos::View<double*> values,
	Connectivity connectivity, IndexType gridIndex)
{
	auto rate = getAppliedRate(gridIndex);
	Kokkos::atomic_sub(&values(connectivity(_reactant, _reactant)), rate);
}

template <typename TNetwork, typename TDerived>
KOKKOS_INLINE_FUNCTION
double
BurstingReaction<TNetwork, TDerived>::computeLeftSideRate(
	ConcentrationsView concentrations, IndexType clusterId, IndexType gridIndex)
{
	return 0.0;
}
} // namespace network
} // namespace core
} // namespace xolotl
