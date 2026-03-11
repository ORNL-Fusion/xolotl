#pragma once

#include <xolotl/util/MathUtils.h>

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
	IndexType reactionId, IndexType cluster0, IndexType cluster1,
	IndexType cluster2) :
	Superclass(reactionData, clusterData, reactionId),
	_reactants({cluster0, cluster1}),
	_product(cluster2)
{
	this->initialize();
}

template <typename TNetwork, typename TDerived>
KOKKOS_INLINE_FUNCTION
NucleationReaction<TNetwork, TDerived>::NucleationReaction(
	ReactionDataRef reactionData, const ClusterData& clusterData,
	IndexType reactionId, const detail::ClusterSet& clusterSet) :
	NucleationReaction(reactionData, clusterData, reactionId,
		clusterSet.cluster0, clusterSet.cluster1, clusterSet.cluster2)
{
}

template <typename TNetwork, typename TDerived>
KOKKOS_INLINE_FUNCTION
double
NucleationReaction<TNetwork, TDerived>::computeRate(
	IndexType gridIndex, double time)
{
	// NE: We say there are 25 bubbles created per fission fragments and there
	// are 2 fission fragments per fission
	// 800H: 71.1 PKA per ion, 90.3 FP produced by cascades above 10 keV per ion
	double rate = (90.3 / 71.1) * this->_clusterData->fissionRate();

	return rate;
}

template <typename TNetwork, typename TDerived>
KOKKOS_INLINE_FUNCTION
void
NucleationReaction<TNetwork, TDerived>::computeConnectivity(
	const Connectivity& connectivity)
{
	// The reactants connects with the reactants
	this->addConnectivity(_reactants[0], _reactants[0], connectivity);
	this->addConnectivity(_reactants[0], _reactants[1], connectivity);
	this->addConnectivity(_reactants[1], _reactants[0], connectivity);
	this->addConnectivity(_reactants[1], _reactants[1], connectivity);
	// The product connects with the reactants
	this->addConnectivity(_product, _reactants[0], connectivity);
	this->addConnectivity(_product, _reactants[1], connectivity);
}

template <typename TNetwork, typename TDerived>
KOKKOS_INLINE_FUNCTION
void
NucleationReaction<TNetwork, TDerived>::computeReducedConnectivity(
	const Connectivity& connectivity)
{
	// The reactants connects with the reactants
	this->addConnectivity(_reactants[0], _reactants[0], connectivity);
	this->addConnectivity(_reactants[1], _reactants[1], connectivity);
	if (_reactants[0] == _reactants[1]) {
		this->addConnectivity(_reactants[0], _reactants[1], connectivity);
		this->addConnectivity(_reactants[1], _reactants[0], connectivity);
	}
	// The product connects with the reactants
	if (_product == _reactants[0])
		this->addConnectivity(_product, _reactants[0], connectivity);
	if (_product == _reactants[1])
		this->addConnectivity(_product, _reactants[1], connectivity);
}

template <typename TNetwork, typename TDerived>
KOKKOS_INLINE_FUNCTION
void
NucleationReaction<TNetwork, TDerived>::computeFlux(
	ConcentrationsView concentrations, FluxesView fluxes, IndexType gridIndex)
{
	// Get the reactant concentrations to know in which regime we are
	double reactConcA = concentrations(_reactants[0]);
	double reactConcB = concentrations(_reactants[1]);

	// Update the concentrations
	if (reactConcA > this->_rate(gridIndex) and
		reactConcB > this->_rate(gridIndex)) {
		Kokkos::atomic_sub(&fluxes(_reactants[0]), this->_rate(gridIndex));
		Kokkos::atomic_sub(&fluxes(_reactants[1]), this->_rate(gridIndex));
		Kokkos::atomic_add(&fluxes(_product), this->_rate(gridIndex));
	}
	else {
		// Get the smallest one
		double minConc = util::min(reactConcA, reactConcB);
		Kokkos::atomic_sub(&fluxes(_reactants[0]), minConc);
		Kokkos::atomic_sub(&fluxes(_reactants[1]), minConc);
		Kokkos::atomic_add(&fluxes(_product), minConc);
	}
}

template <typename TNetwork, typename TDerived>
KOKKOS_INLINE_FUNCTION
void
NucleationReaction<TNetwork, TDerived>::computePartialDerivatives(
	ConcentrationsView concentrations, Kokkos::View<double*> values,
	IndexType gridIndex)
{
	// Get the reactant concentrations to know in which regime we are
	double reactConcA = concentrations(_reactants[0]);
	double reactConcB = concentrations(_reactants[1]);

	// Update the concentrations
	if (reactConcA > this->_rate(gridIndex) and
		reactConcB > this->_rate(gridIndex)) {
		// Nothing
	}
	else {
		// Need to know which one is smallest
		if (reactConcA < reactConcB) {
			Kokkos::atomic_sub(&values(_connEntries[0][0]), 1.0);
			Kokkos::atomic_sub(&values(_connEntries[1][0]), 1.0);
			Kokkos::atomic_add(&values(_connEntries[2][0]), 1.0);
		}
		else {
			Kokkos::atomic_sub(&values(_connEntries[0][1]), 1.0);
			Kokkos::atomic_sub(&values(_connEntries[1][1]), 1.0);
			Kokkos::atomic_add(&values(_connEntries[2][1]), 1.0);
		}
	}
}

template <typename TNetwork, typename TDerived>
KOKKOS_INLINE_FUNCTION
void
NucleationReaction<TNetwork, TDerived>::computeReducedPartialDerivatives(
	ConcentrationsView concentrations, Kokkos::View<double*> values,
	IndexType gridIndex)
{
	// Get the reactant concentrations to know in which regime we are
	double reactConcA = concentrations(_reactants[0]);
	double reactConcB = concentrations(_reactants[1]);

	// Update the concentrations
	if (reactConcA > this->_rate(gridIndex) and
		reactConcB > this->_rate(gridIndex)) {
		// Nothing
	}
	else {
		// Need to know which one is smallest
		if (reactConcA < reactConcB) {
			Kokkos::atomic_sub(&values(_connEntries[0][0]), 1.0);
			if (_reactants[0] == _reactants[1])
				Kokkos::atomic_sub(&values(_connEntries[1][0]), 1.0);
			if (_reactants[0] == _product)
				Kokkos::atomic_add(&values(_connEntries[2][0]), 1.0);
		}
		else {
			if (_reactants[0] == _reactants[1])
				Kokkos::atomic_sub(&values(_connEntries[0][1]), 1.0);
			Kokkos::atomic_sub(&values(_connEntries[1][1]), 1.0);
			if (_reactants[1] == _product)
				Kokkos::atomic_add(&values(_connEntries[2][1]), 1.0);
		}
	}
}

template <typename TNetwork, typename TDerived>
KOKKOS_INLINE_FUNCTION
void
NucleationReaction<TNetwork, TDerived>::mapJacobianEntries(
	Connectivity connectivity)
{
	_connEntries[0][0] = connectivity(_reactants[0], _reactants[0]);
	_connEntries[0][1] = connectivity(_reactants[0], _reactants[1]);
	_connEntries[1][0] = connectivity(_reactants[1], _reactants[0]);
	_connEntries[1][1] = connectivity(_reactants[1], _reactants[1]);
	_connEntries[2][0] = connectivity(_product, _reactants[0]);
	_connEntries[2][1] = connectivity(_product, _reactants[1]);
}
} // namespace network
} // namespace core
} // namespace xolotl
