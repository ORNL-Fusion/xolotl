#pragma once

#include <xolotl/core/network/Reaction.h>
#include <xolotl/core/network/SpeciesEnumSequence.h>

namespace xolotl
{
namespace core
{
namespace network
{
template <typename TNetwork, typename TDerived>
class NucleationReaction : public Reaction<TNetwork, TDerived>
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
	using AmountType = typename Superclass::AmountType;
	using ReactionDataRef = typename Superclass::ReactionDataRef;

	NucleationReaction() = default;

	KOKKOS_INLINE_FUNCTION
	NucleationReaction(ReactionDataRef reactionData, ClusterDataRef clusterData,
		IndexType reactionId, IndexType cluster0, IndexType cluster1);

	KOKKOS_INLINE_FUNCTION
	NucleationReaction(ReactionDataRef reactionData, ClusterDataRef clusterData,
		IndexType reactionId, const detail::ClusterSet& clusterSet);

	static detail::CoefficientsView allocateCoefficientsView(IndexType)
	{
		return detail::CoefficientsView();
	}

private:
	KOKKOS_INLINE_FUNCTION
	void
	computeCoefficients()
	{
		// No coefs
	}

	KOKKOS_INLINE_FUNCTION
	double
	computeRate(IndexType gridIndex);

	KOKKOS_INLINE_FUNCTION
	void
	computeConnectivity(const Connectivity& connectivity);

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
	double
	computeLeftSideRate(ConcentrationsView concentrations, IndexType clusterId,
		IndexType gridIndex)
	{
		return 0.0;
	}

protected:
	IndexType _reactant;
	IndexType _product;
	static constexpr auto invalidIndex = Superclass::invalidIndex;
};
} // namespace network
} // namespace core
} // namespace xolotl

#include <xolotl/core/network/detail/NucleationReactionGenerator.h>
