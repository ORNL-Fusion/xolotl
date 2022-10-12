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
	computeConnectivity(const Connectivity& connectivity)
	{
		// Nothing
	}

	KOKKOS_INLINE_FUNCTION
	void
	computeReducedConnectivity(const Connectivity& connectivity)
	{
		// Nothing
	}

	KOKKOS_INLINE_FUNCTION
	void
	computeFlux(ConcentrationsView concentrations, FluxesView fluxes,
		IndexType gridIndex)
	{
		Kokkos::atomic_sub(&fluxes(_reactant), this->_rate(gridIndex));
		Kokkos::atomic_add(&fluxes(_product), this->_rate(gridIndex));
	}

	KOKKOS_INLINE_FUNCTION
	void
	computePartialDerivatives(ConcentrationsView concentrations,
		Kokkos::View<double*> values, IndexType gridIndex)
	{
		return;
	}

	KOKKOS_INLINE_FUNCTION
	void
	computeReducedPartialDerivatives(ConcentrationsView concentrations,
		Kokkos::View<double*> values, IndexType gridIndex)
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
		// Nothing
	}

protected:
	IndexType _reactant;
	IndexType _product;
	static constexpr auto invalidIndex = Superclass::invalidIndex;
};
} // namespace network
} // namespace core
} // namespace xolotl

#include <xolotl/core/network/detail/TransformReactionGenerator.h>
