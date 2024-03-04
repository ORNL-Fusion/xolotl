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
 * @brief Class implementing bursting reaction where
 * He_mV_n -> V_n with a given rate
 * and He_m -> 0
 *
 * @tparam TNetwork The network type
 * @tparam TDerived The derived class type.
 */
template <typename TNetwork, typename TDerived>
class BurstingReaction : public Reaction<TNetwork, TDerived>
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
	using Composition = typename Superclass::Composition;
	using Region = typename Superclass::Region;
	using ConnectivitiesView = typename Superclass::ConnectivitiesView;
	using ConnectivitiesPairView = typename Superclass::ConnectivitiesPairView;
	using BelongingView = typename Superclass::BelongingView;
	using OwnedSubMapView = typename Superclass::OwnedSubMapView;
	using AmountType = typename Superclass::AmountType;
	using ReactionDataRef = typename Superclass::ReactionDataRef;
	using ClusterData = typename Superclass::ClusterData;

	BurstingReaction() = default;

	KOKKOS_INLINE_FUNCTION
	BurstingReaction(ReactionDataRef reactionData,
		const ClusterData& clusterData, IndexType reactionId,
		IndexType cluster0, IndexType cluster1);

	KOKKOS_INLINE_FUNCTION
	BurstingReaction(ReactionDataRef reactionData,
		const ClusterData& clusterData, IndexType reactionId,
		const detail::ClusterSet& clusterSet);

	static detail::CoefficientsView
	allocateCoefficientsView(IndexType size)
	{
		return detail::CoefficientsView("Bursting Coefficients", size,
			Superclass::coeffsSingleExtent, 1, 1,
			Superclass::coeffsSingleExtent);
	}

	static detail::ConstantRateView
	allocateConstantRateView(IndexType, IndexType)
	{
		return detail::ConstantRateView();
	}

	using Superclass::updateRates;

	KOKKOS_INLINE_FUNCTION
	void
	updateLargestRates(double largestRate);

	KOKKOS_INLINE_FUNCTION
	void
	computeCoefficients();

	KOKKOS_INLINE_FUNCTION
	void
	computeFlux(ConcentrationsView concentrations, FluxesView fluxes,
		IndexType gridIndex);

	KOKKOS_INLINE_FUNCTION
	void
	computePartialDerivatives(ConcentrationsView concentrations,
		Kokkos::View<double*> values, IndexType gridIndex);

private:
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
	computeReducedPartialDerivatives(ConcentrationsView concentrations,
		Kokkos::View<double*> values, IndexType gridIndex);

	KOKKOS_INLINE_FUNCTION
	void
	computeConstantRates(ConcentrationsView concentrations, RatesView rates,
		BelongingView isInSub, IndexType subId, IndexType gridIndex)
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
		IndexType gridIndex);

	KOKKOS_INLINE_FUNCTION
	void
	mapJacobianEntries(Connectivity connectivity);

	KOKKOS_INLINE_FUNCTION
	void
	mapRateEntries(ConnectivitiesPairView connectivityRow,
		ConnectivitiesPairView connectivityEntries, BelongingView isInSub,
		OwnedSubMapView backMap, IndexType subId)
	{
		return;
	}

protected:
	IndexType _reactant;
	double _reactantVolume;
	IndexType _product;
	static constexpr auto invalidIndex = Superclass::invalidIndex;

	static constexpr auto nMomentIds = Superclass::nMomentIds;
	util::Array<IndexType, nMomentIds> _reactantMomentIds;

	util::Array<IndexType, 2, 1 + nMomentIds, 1, 1 + nMomentIds> _connEntries;
};
} // namespace network
} // namespace core
} // namespace xolotl

#include <xolotl/core/network/detail/BurstingReactionGenerator.h>
