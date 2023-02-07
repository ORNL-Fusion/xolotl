#pragma once

#include <xolotl/core/network/Reaction.h>

namespace xolotl
{
namespace core
{
namespace network
{
/**
 * @brief Class implementing re-solution reaction where
 * X_n -> X_(n-1) + X_1 with a given rate
 *
 * @tparam TNetwork The network type
 * @tparam TDerived The derived class type.
 */
template <typename TNetwork, typename TDerived>
class PartialReSolutionReaction : public Reaction<TNetwork, TDerived>
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
	using Composition = typename Superclass::Composition;
	using Region = typename Superclass::Region;
	using AmountType = typename Superclass::AmountType;
	using ReactionDataRef = typename Superclass::ReactionDataRef;
	using ClusterData = typename Superclass::ClusterData;
	using ReflectedRegion = typename Superclass::ReflectedRegion;

	PartialReSolutionReaction() = default;

	KOKKOS_INLINE_FUNCTION
	PartialReSolutionReaction(ReactionDataRef reactionData,
		const ClusterData& clusterData, IndexType reactionId,
		IndexType cluster0, IndexType cluster1, IndexType cluster2);

	KOKKOS_INLINE_FUNCTION
	PartialReSolutionReaction(ReactionDataRef reactionData,
		const ClusterData& clusterData, IndexType reactionId,
		const detail::ClusterSet& clusterSet);

	static detail::CoefficientsView
	allocateCoefficientsView(IndexType size)
	{
		return detail::CoefficientsView("PartialReSolution Coefficients", size,
			Superclass::coeffsSingleExtent, 1, 3,
			Superclass::coeffsSingleExtent);
	}

private:
	KOKKOS_INLINE_FUNCTION
	void
	computeCoefficients();

	KOKKOS_INLINE_FUNCTION
	double
	computeRate(IndexType gridIndex, double time = 0.0);

	KOKKOS_INLINE_FUNCTION
	void
	computeConnectivity(const Connectivity& connectivity);

	KOKKOS_INLINE_FUNCTION
	void
	computeReducedConnectivity(const Connectivity& connectivity);

	KOKKOS_INLINE_FUNCTION
	void
	computeFlux(ConcentrationsView concentrations, FluxesView fluxes,
		IndexType gridIndex);

	KOKKOS_INLINE_FUNCTION
	void
	computePartialDerivatives(ConcentrationsView concentrations,
		Kokkos::View<double*> values, IndexType gridIndex);

	KOKKOS_INLINE_FUNCTION
	void
	computeReducedPartialDerivatives(ConcentrationsView concentrations,
		Kokkos::View<double*> values, IndexType gridIndex);

	KOKKOS_INLINE_FUNCTION
	void
	computeConstantRates(ConcentrationsView concentrations, RatesView rates,
		BelongingView isInSub, OwnedSubMapView backMap, IndexType gridIndex);

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
		IndexType gridIndex);

	KOKKOS_INLINE_FUNCTION
	void
	mapJacobianEntries(Connectivity connectivity);

protected:
	IndexType _reactant;
	double _reactantVolume;
	static constexpr auto invalidIndex = Superclass::invalidIndex;
	util::Array<IndexType, 2> _products{invalidIndex, invalidIndex};
	util::Array<double, 2> _productVolumes{0.0, 0.0};

	static constexpr auto nMomentIds = Superclass::nMomentIds;
	util::Array<IndexType, nMomentIds> _reactantMomentIds;
	util::Array<IndexType, 2, nMomentIds> _productMomentIds;

	util::Array<IndexType, 3, 1 + nMomentIds, 1, 1 + nMomentIds> _connEntries;
};
} // namespace network
} // namespace core
} // namespace xolotl

#include <xolotl/core/network/detail/PartialReSolutionReactionGenerator.h>
