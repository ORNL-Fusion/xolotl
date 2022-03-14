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
 * @brief Class implementing reaction where the other cluster is missing from
 * the network but its rate is accounted for as a constant rate that can be
 * updated.
 *
 * @tparam TNetwork The network type
 * @tparam TDerived The derived class type.
 */
template <typename TNetwork, typename TDerived>
class ConstantReaction : public Reaction<TNetwork, TDerived>
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
	using BelongingView = typename Superclass::BelongingView;
	using OwnedSubMapView = typename Superclass::OwnedSubMapView;
	using AmountType = typename Superclass::AmountType;
	using ReactionDataRef = typename Superclass::ReactionDataRef;
	using ClusterData = typename Superclass::ClusterData;
	using RateVector = IReactionNetwork::RateVector;

	ConstantReaction() = default;

	KOKKOS_INLINE_FUNCTION
	ConstantReaction(ReactionDataRef reactionData,
		const ClusterData& clusterData, IndexType reactionId,
		IndexType cluster0, IndexType cluster1) :
		Superclass(reactionData, clusterData, reactionId),
		_reactants({cluster0, cluster1})
	{
		this->initialize();
	}

	KOKKOS_INLINE_FUNCTION
	ConstantReaction(ReactionDataRef reactionData,
		const ClusterData& clusterData, IndexType reactionId,
		const detail::ClusterSet& clusterSet) :
		ConstantReaction(reactionData, clusterData, reactionId,
			clusterSet.cluster0, clusterSet.cluster1)
	{
	}

	static detail::CoefficientsView allocateCoefficientsView(IndexType)
	{
		return detail::CoefficientsView();
	}

	KOKKOS_INLINE_FUNCTION
	double
	computeRate(IndexType gridIndex)
	{
		return _constantRates[0][0][0][0];
	}

	KOKKOS_INLINE_FUNCTION
	void
	setRate(RatesView rates)
	{
		auto dof = rates.extent(0);
		if (_reactants[1] == invalidIndex)
			_constantRates[0][0][0][0] = rates(_reactants[0], dof);
		else
			_constantRates[0][0][0][0] = rates(_reactants[0], _reactants[1]);
	}

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
		if (_reactants[1] == invalidIndex)
			return;
		// The second reactant connects with the first reactant
		this->addConnectivity(_reactants[0], _reactants[1], connectivity);
	}

	KOKKOS_INLINE_FUNCTION
	void
	computeReducedConnectivity(const Connectivity& connectivity)
	{
		// The second reactant connects with the first reactant
		if (_reactants[0] == _reactants[1])
			this->addConnectivity(_reactants[0], _reactants[1], connectivity);
	}

	KOKKOS_INLINE_FUNCTION
	void
	computeFlux(ConcentrationsView concentrations, FluxesView fluxes,
		IndexType gridIndex)
	{
		if (_reactants[1] == invalidIndex)
			Kokkos::atomic_add(&fluxes(_reactants[0]), this->_rate(gridIndex));
		else
			Kokkos::atomic_add(&fluxes(_reactants[0]),
				concentrations(_reactants[1]) * this->_rate(gridIndex));
	}

	KOKKOS_INLINE_FUNCTION
	void
	computePartialDerivatives(ConcentrationsView concentrations,
		Kokkos::View<double*> values, IndexType gridIndex)
	{
		if (_reactants[1] == invalidIndex)
			return;
		Kokkos::atomic_add(
			&values(_connEntries[0][0][0][0]), this->_rate(gridIndex));
	}

	KOKKOS_INLINE_FUNCTION
	void
	computeReducedPartialDerivatives(ConcentrationsView concentrations,
		Kokkos::View<double*> values, IndexType gridIndex)
	{
		if (_reactants[0] == _reactants[1])
			Kokkos::atomic_add(
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
		if (_reactants[1] == invalidIndex)
			return;
		_connEntries[0][0][0][0] = connectivity(_reactants[0], _reactants[1]);
	}

protected:
	static constexpr auto invalidIndex = Superclass::invalidIndex;
	util::Array<IndexType, 2> _reactants{invalidIndex, invalidIndex};

	util::Array<IndexType, 1, 1, 1, 1> _connEntries;
	util::Array<double, 1, 1, 1, 1> _constantRates;
};
} // namespace network
} // namespace core
} // namespace xolotl

#include <xolotl/core/network/detail/ConstantReactionGenerator.h>
