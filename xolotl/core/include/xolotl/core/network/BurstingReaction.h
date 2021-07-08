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
	using ClusterDataRef = typename Superclass::ClusterDataRef;
	using IndexType = typename Superclass::IndexType;
	using Connectivity = typename Superclass::Connectivity;
	using ConcentrationsView = typename Superclass::ConcentrationsView;
	using FluxesView = typename Superclass::FluxesView;
	using Composition = typename Superclass::Composition;
	using Region = typename Superclass::Region;
	using AmountType = typename Superclass::AmountType;
	using ReactionDataRef = typename Superclass::ReactionDataRef;

	BurstingReaction() = default;

	KOKKOS_INLINE_FUNCTION
	BurstingReaction(ReactionDataRef reactionData, ClusterDataRef clusterData,
		IndexType reactionId, IndexType cluster0, IndexType cluster1);

	KOKKOS_INLINE_FUNCTION
	BurstingReaction(ReactionDataRef reactionData, ClusterDataRef clusterData,
		IndexType reactionId, const detail::ClusterSet& clusterSet);

	static detail::CoefficientsView
	allocateCoefficientsView(IndexType size)
	{
		return detail::CoefficientsView("Bursting Coefficients", size,
			Superclass::coeffsSingleExtent, 1, 1,
			Superclass::coeffsSingleExtent);
	}

	using Superclass::updateRates;

	KOKKOS_INLINE_FUNCTION
	void
	updateRates(double largestRate);

private:
	KOKKOS_INLINE_FUNCTION
	void
	computeCoefficients();

	KOKKOS_INLINE_FUNCTION
	double
	computeRate(IndexType gridIndex);

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
	void
	computeFlux(ConcentrationsView concentrations, FluxesView fluxes,
		IndexType gridIndex);

	KOKKOS_INLINE_FUNCTION
	void
	computePartialDerivatives(ConcentrationsView concentrations,
		Kokkos::View<double*> values, Connectivity connectivity,
		IndexType gridIndex);

	KOKKOS_INLINE_FUNCTION
	void
	computeReducedPartialDerivatives(ConcentrationsView concentrations,
		Kokkos::View<double*> values, Connectivity connectivity,
		IndexType gridIndex);

	KOKKOS_INLINE_FUNCTION
	double
	computeLeftSideRate(ConcentrationsView concentrations, IndexType clusterId,
		IndexType gridIndex);

protected:
	IndexType _reactant;
	IndexType _product;
	static constexpr auto invalidIndex = Superclass::invalidIndex;

	static constexpr auto nMomentIds = Superclass::nMomentIds;
	Kokkos::Array<IndexType, nMomentIds> _reactantMomentIds;
};
} // namespace network
} // namespace core
} // namespace xolotl

#include <xolotl/core/network/detail/BurstingReactionGenerator.h>
