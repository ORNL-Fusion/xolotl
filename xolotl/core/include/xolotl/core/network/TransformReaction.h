#pragma once

#include <xolotl/core/network/Reaction.h>
#include <xolotl/core/network/SpeciesEnumSequence.h>

namespace xolotl
{
namespace core
{
namespace network
{
/**
 * @brief Class implementing transform reaction where
 * reactant -> product with a given constant rate
 *
 * @tparam TNetwork The network type
 * @tparam TDerived The derived class type.
 */
template <typename TNetwork, typename TDerived>
class TransformReaction : public Reaction<TNetwork, TDerived>
{
	friend class Reaction<TNetwork, TDerived>;

public:
	using NetworkType = TNetwork;
	using Superclass = Reaction<TNetwork, TDerived>;
	using IndexType = typename Superclass::IndexType;
	using Connectivity = typename Superclass::Connectivity;
	using ConcentrationsView = typename Superclass::ConcentrationsView;
	using FluxesView = typename Superclass::FluxesView;
	using RatesView = typename Superclass::RatesView;
	using ConnectivitiesView = typename Superclass::ConnectivitiesView;
	using BelongingView = typename Superclass::BelongingView;
	using OwnedSubMapView = typename Superclass::OwnedSubMapView;
	using AmountType = typename Superclass::AmountType;
	using ReactionDataRef = typename Superclass::ReactionDataRef;
	using ClusterData = typename Superclass::ClusterData;

	TransformReaction() = default;

	KOKKOS_INLINE_FUNCTION
	TransformReaction(ReactionDataRef reactionData,
		const ClusterData& clusterData, IndexType reactionId,
		IndexType cluster0, IndexType cluster1) :
		Superclass(reactionData, clusterData, reactionId),
		_reactant(cluster0),
		_product(cluster1)
	{
		this->initialize();
	}

	KOKKOS_INLINE_FUNCTION
	TransformReaction(ReactionDataRef reactionData,
		const ClusterData& clusterData, IndexType reactionId,
		const detail::ClusterSet& clusterSet) :
		TransformReaction(reactionData, clusterData, reactionId,
			clusterSet.cluster0, clusterSet.cluster1)
	{
	}

	static detail::CoefficientsView
	allocateCoefficientsView(IndexType)
	{
		return detail::CoefficientsView();
	}

	KOKKOS_INLINE_FUNCTION
	double
	computeRate(IndexType gridIndex, double time = 0.0);

private:
	KOKKOS_INLINE_FUNCTION
	void
	computeCoefficients()
	{
		// No coefs
	}

	KOKKOS_INLINE_FUNCTION
	void
	computeConnectivity(const Connectivity& connectivity)
	{
		// Everything connects with the reactant
		this->addConnectivity(_reactant, _reactant, connectivity);
		this->addConnectivity(_product, _reactant, connectivity);
	}

	KOKKOS_INLINE_FUNCTION
	void
	computeReducedConnectivity(const Connectivity& connectivity)
	{
		// Everything connects with the reactant
		this->addConnectivity(_reactant, _reactant, connectivity);
		if (_product == _reactant)
			this->addConnectivity(_product, _reactant, connectivity);
	}

	KOKKOS_INLINE_FUNCTION
	void
	computeFlux(ConcentrationsView concentrations, FluxesView fluxes,
		IndexType gridIndex)
	{
		Kokkos::atomic_sub(&fluxes(_reactant),
			this->_rate(gridIndex) * concentrations(_reactant));
		Kokkos::atomic_add(&fluxes(_product),
			this->_rate(gridIndex) * concentrations(_reactant));
	}

	KOKKOS_INLINE_FUNCTION
	void
	computePartialDerivatives(ConcentrationsView concentrations,
		Kokkos::View<double*> values, IndexType gridIndex)
	{
		Kokkos::atomic_sub(
			&values(_connEntries[0][0][0][0]), this->_rate(gridIndex));
		Kokkos::atomic_add(
			&values(_connEntries[1][0][0][0]), this->_rate(gridIndex));
	}

	KOKKOS_INLINE_FUNCTION
	void
	computeReducedPartialDerivatives(ConcentrationsView concentrations,
		Kokkos::View<double*> values, IndexType gridIndex)
	{
		Kokkos::atomic_sub(
			&values(_connEntries[0][0][0][0]), this->_rate(gridIndex));
		if (_product == _reactant)
			Kokkos::atomic_add(
				&values(_connEntries[1][0][0][0]), this->_rate(gridIndex));
	}

	KOKKOS_INLINE_FUNCTION
	void
	computeConstantRates(ConcentrationsView concentrations, RatesView rates,
		BelongingView isInSub, OwnedSubMapView backMap, IndexType gridIndex)
	{
		return;
	}

	KOKKOS_INLINE_FUNCTION
	void
	getConstantConnectivities(ConnectivitiesView conns, BelongingView isInSub,
		OwnedSubMapView backMap)
	{
		return;
	}

	KOKKOS_INLINE_FUNCTION
	double
	computeLeftSideRate(ConcentrationsView concentrations, IndexType clusterId,
		IndexType gridIndex)
	{
		return 0.0;
	}

	KOKKOS_INLINE_FUNCTION
	double
	computeNetSigma(ConcentrationsView concentrations, IndexType clusterId,
		IndexType gridIndex)
	{
		return 0.0;
	}

	KOKKOS_INLINE_FUNCTION
	void
	mapJacobianEntries(Connectivity connectivity)
	{
		_connEntries[0][0][0][0] = connectivity(_reactant, _reactant);
		_connEntries[1][0][0][0] = connectivity(_product, _reactant);
	}

protected:
	IndexType _reactant;
	IndexType _product;
	static constexpr auto invalidIndex = Superclass::invalidIndex;

	util::Array<IndexType, 2, 1, 1, 1> _connEntries;
};
} // namespace network
} // namespace core
} // namespace xolotl

#include <xolotl/core/network/detail/TransformReactionGenerator.h>
