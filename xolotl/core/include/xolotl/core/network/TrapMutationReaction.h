#pragma once

#include <xolotl/core/network/Reaction.h>

namespace xolotl
{
namespace core
{
namespace network
{
template <typename TNetwork, typename TDerived>
class TrapMutationReaction : public Reaction<TNetwork, TDerived>
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

protected:
	static constexpr auto invalidIndex = Superclass::invalidIndex;

public:
	TrapMutationReaction() = default;

	KOKKOS_INLINE_FUNCTION
	TrapMutationReaction(ReactionDataRef reactionData,
		ClusterDataRef clusterData, IndexType reactionId, IndexType cluster0,
		IndexType cluster1, IndexType cluster2);

	KOKKOS_INLINE_FUNCTION
	TrapMutationReaction(ReactionDataRef reactionData,
		ClusterDataRef clusterData, IndexType reactionId,
		const detail::ClusterSet& clusterSet);

	static detail::CoefficientsView allocateCoefficientsView(IndexType)
	{
		return detail::CoefficientsView();
	}

	KOKKOS_INLINE_FUNCTION
	double
	computeRate(double largestRate);

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

private:
	AmountType _heAmount;
	IndexType _heClId;
	IndexType _heVClId;
	IndexType _iClId;
};
} // namespace network
} // namespace core
} // namespace xolotl

#include <xolotl/core/network/detail/TrapMutationReactionGenerator.h>
