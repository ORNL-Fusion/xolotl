#pragma once

#include <math.h>

#include <xolotl/core/Constants.h>
#include <xolotl/core/network/Reaction.h>
#include <xolotl/core/network/SpeciesEnumSequence.h>

namespace xolotl
{
namespace core
{
namespace network
{
/**
 * @brief Class implementing sink reaction where
 * reactant -> 0 with a given rate
 *
 * @tparam TNetwork The network type
 * @tparam TDerived The derived class type.
 */
template <typename TNetwork, typename TDerived>
class SinkReaction : public Reaction<TNetwork, TDerived>
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

	SinkReaction() = default;

	KOKKOS_INLINE_FUNCTION
	SinkReaction(ReactionDataRef reactionData, const ClusterData& clusterData,
		IndexType reactionId, IndexType cluster0) :
		Superclass(reactionData, clusterData, reactionId),
		_reactant(cluster0)
	{
		this->initialize();
	}

	KOKKOS_INLINE_FUNCTION
	SinkReaction(ReactionDataRef reactionData, const ClusterData& clusterData,
		IndexType reactionId, const detail::ClusterSet& clusterSet) :
		SinkReaction(reactionData, clusterData, reactionId, clusterSet.cluster0)
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
		// The reactant connects with the reactant
		this->addConnectivity(_reactant, _reactant, connectivity);
	}

	KOKKOS_INLINE_FUNCTION
	void
	computeReducedConnectivity(const Connectivity& connectivity)
	{
		// The reactant connects with the reactant
		this->addConnectivity(_reactant, _reactant, connectivity);
	}

	KOKKOS_INLINE_FUNCTION
	void
	computeFlux(ConcentrationsView concentrations, FluxesView fluxes,
		IndexType gridIndex)
	{
		Kokkos::atomic_sub(&fluxes(_reactant),
			concentrations(_reactant) * this->_rate(gridIndex));
	}

	KOKKOS_INLINE_FUNCTION
	void
	computePartialDerivatives(ConcentrationsView concentrations,
		Kokkos::View<double*> values, IndexType gridIndex)
	{
		Kokkos::atomic_sub(
			&values(_connEntries[0][0][0][0]), this->_rate(gridIndex));
	}

	KOKKOS_INLINE_FUNCTION
	void
	computeReducedPartialDerivatives(ConcentrationsView concentrations,
		Kokkos::View<double*> values, IndexType gridIndex)
	{
		Kokkos::atomic_sub(
			&values(_connEntries[0][0][0][0]), this->_rate(gridIndex));
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
	void
	mapJacobianEntries(Connectivity connectivity)
	{
		_connEntries[0][0][0][0] = connectivity(_reactant, _reactant);
	}

protected:
	IndexType _reactant;
	static constexpr auto invalidIndex = Superclass::invalidIndex;

	util::Array<IndexType, 1, 1, 1, 1> _connEntries;
};
} // namespace network
} // namespace core
} // namespace xolotl

#include <xolotl/core/network/detail/SinkReactionGenerator.h>
