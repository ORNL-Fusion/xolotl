#pragma once

#include <xolotl/core/network/Reaction.h>

namespace xolotl
{
namespace core
{
namespace network
{
/**
 * @brief Class implementing near surface trap mutation reactions where
 * He_n -> He_nV_m + I_m
 * with a given rate at specific depths.
 *
 * The depths depend on the material orientation and temperature.
 *
 * A cluster can also desorb, meaning that the rate of the trap mutation
 * reaction will depend on the "left side" rate.
 *
 * If the attenuation is ON, the rates will decrease as a function of the
 * concentration value given to
 * IPSIReactionNetwork::updateTrapMutationDisappearingRate(conc).
 *
 * @tparam TNetwork The network type
 * @tparam TDerived The derived class type.
 */
template <typename TNetwork, typename TDerived>
class TrapMutationReaction : public Reaction<TNetwork, TDerived>
{
	friend class Reaction<TNetwork, TDerived>;

public:
	using NetworkType = TNetwork;
	using Superclass = Reaction<TNetwork, TDerived>;
	using Species = typename Superclass::Species;
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

protected:
	static constexpr auto invalidIndex = Superclass::invalidIndex;

public:
	TrapMutationReaction() = default;

	KOKKOS_INLINE_FUNCTION
	TrapMutationReaction(ReactionDataRef reactionData,
		const ClusterData& clusterData, IndexType reactionId,
		IndexType cluster0, IndexType cluster1, IndexType cluster2);

	KOKKOS_INLINE_FUNCTION
	TrapMutationReaction(ReactionDataRef reactionData,
		const ClusterData& clusterData, IndexType reactionId,
		const detail::ClusterSet& clusterSet);

	static detail::CoefficientsView
	allocateCoefficientsView(IndexType)
	{
		return detail::CoefficientsView();
	}

	static detail::ConstantRateView allocateConstantRateView(
		IndexType, IndexType)
	{
		return detail::ConstantRateView();
	}

	KOKKOS_INLINE_FUNCTION
	double
	computeRate(double largestRate);

	using Superclass::updateRates;

	KOKKOS_INLINE_FUNCTION
	void
	updateRates(double largestRate);

private:
	KOKKOS_INLINE_FUNCTION
	void
	computeCoefficients()
	{
		// No coefs
	}

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
	double
	getAppliedRate(IndexType gridIndex) const;

	KOKKOS_INLINE_FUNCTION
	bool
	getEnabled() const;

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
	mapJacobianEntries(Connectivity connectivity)
	{
		_connEntries[0][0][0][0] = connectivity(_heClId, _heClId);
		_connEntries[1][0][0][0] = connectivity(_heVClId, _heClId);
		_connEntries[2][0][0][0] = connectivity(_iClId, _heClId);
	}

private:
	IndexType _heClId;
	IndexType _heVClId;
	IndexType _iClId;
	AmountType _heAmount{};
	AmountType _vSize{};

	util::Array<IndexType, 3, 1, 1, 1> _connEntries;
};
} // namespace network
} // namespace core
} // namespace xolotl

#include <xolotl/core/network/detail/TrapMutationReactionGenerator.h>
