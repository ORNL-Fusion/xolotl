#pragma once

namespace xolotl
{
namespace core
{
namespace network
{
template <typename TNetwork, typename TDerived>
KOKKOS_INLINE_FUNCTION
NucleationReaction<TNetwork, TDerived>::NucleationReaction(
	ReactionDataRef reactionData, const ClusterData& clusterData,
	IndexType reactionId, IndexType cluster0, IndexType cluster1) :
	Superclass(reactionData, clusterData, reactionId),
	_reactant(cluster0),
	_product(cluster1)
{
	this->initialize();
}

template <typename TNetwork, typename TDerived>
KOKKOS_INLINE_FUNCTION
NucleationReaction<TNetwork, TDerived>::NucleationReaction(
	ReactionDataRef reactionData, const ClusterData& clusterData,
	IndexType reactionId, const detail::ClusterSet& clusterSet) :
	NucleationReaction(reactionData, clusterData, reactionId,
		clusterSet.cluster0, clusterSet.cluster1)
{
}

template <typename TNetwork, typename TDerived>
KOKKOS_INLINE_FUNCTION
double
NucleationReaction<TNetwork, TDerived>::computeRate(
	IndexType gridIndex, double time)
{
	// We say there are 25 bubbles created per fission fragments and there
	// are 2 fission fragments per fission
	double rate = 50.0 * this->_clusterData->fissionRate();

	return rate;
}

template <typename TNetwork, typename TDerived>
KOKKOS_INLINE_FUNCTION
void
NucleationReaction<TNetwork, TDerived>::computeConnectivity(
	const Connectivity& connectivity)
{
	// The reactant connects with the reactant
	this->addConnectivity(_reactant, _reactant, connectivity);
	// The product connects with the reactant
	this->addConnectivity(_product, _reactant, connectivity);
}

template <typename TNetwork, typename TDerived>
KOKKOS_INLINE_FUNCTION
void
NucleationReaction<TNetwork, TDerived>::computeReducedConnectivity(
	const Connectivity& connectivity)
{
	// The reactant connects with the reactant
	this->addConnectivity(_reactant, _reactant, connectivity);
	// The product connects with the reactant
	if (_product == _reactant)
		this->addConnectivity(_product, _reactant, connectivity);
}

template <typename TNetwork, typename TDerived>
KOKKOS_INLINE_FUNCTION
void
NucleationReaction<TNetwork, TDerived>::computeFlux(
	ConcentrationsView concentrations, FluxesView fluxes, IndexType gridIndex)
{
	// Get the single concentration to know in which regime we are
	double singleConc = concentrations(_reactant);

	// Update the concentrations
	if (singleConc > 2.0 * this->_rate(gridIndex)) {
		Kokkos::atomic_sub(&fluxes(_reactant), 2.0 * this->_rate(gridIndex));
		Kokkos::atomic_add(&fluxes(_product), this->_rate(gridIndex));
	}
	else {
		Kokkos::atomic_sub(&fluxes(_reactant), singleConc);
		Kokkos::atomic_add(&fluxes(_product), singleConc / 2.0);
	}
}

template <typename TNetwork, typename TDerived>
KOKKOS_INLINE_FUNCTION
void
NucleationReaction<TNetwork, TDerived>::computePartialDerivatives(
	ConcentrationsView concentrations, Kokkos::View<double*> values,
	IndexType gridIndex)
{
	// Get the single concentration to know in which regime we are
	double singleConc = concentrations(_reactant);

	// Update the partials
	if (singleConc > 2.0 * this->_rate(gridIndex)) {
		// Nothing
	}
	else {
		Kokkos::atomic_sub(&values(_connEntries[0][0][0][0]), 1.0);
		Kokkos::atomic_add(&values(_connEntries[1][0][0][0]), 0.5);
	}
}

template <typename TNetwork, typename TDerived>
KOKKOS_INLINE_FUNCTION
void
NucleationReaction<TNetwork, TDerived>::computeReducedPartialDerivatives(
	ConcentrationsView concentrations, Kokkos::View<double*> values,
	IndexType gridIndex)
{
	// Get the single concentration to know in which regime we are
	double singleConc = concentrations(_reactant);

	// Update the partials
	if (singleConc > 2.0 * this->_rate(gridIndex)) {
		// Nothing
	}
	else {
		Kokkos::atomic_sub(&values(_connEntries[0][0][0][0]), 1.0);
		if (_product == _reactant)
			Kokkos::atomic_add(&values(_connEntries[1][0][0][0]), 0.5);
	}
}

template <typename TNetwork, typename TDerived>
KOKKOS_INLINE_FUNCTION
void
NucleationReaction<TNetwork, TDerived>::mapJacobianEntries(
	Connectivity connectivity)
{
	_connEntries[0][0][0][0] = connectivity(_reactant, _reactant);
	_connEntries[1][0][0][0] = connectivity(_product, _reactant);
}
} // namespace network
} // namespace core
} // namespace xolotl
