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
 * @brief Class implementing Trap reaction where
 * reactant_mobile <-> reactant_trapped with a given rate
 * of trapping and detrapping
 *
 * @tparam TNetwork The network type
 * @tparam TDerived The derived class type.
 */
template <typename TNetwork, typename TDerived>
class TrapReaction : public Reaction<TNetwork, TDerived>
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
	using ConnectivitiesPairView = typename Superclass::ConnectivitiesPairView;
	using BelongingView = typename Superclass::BelongingView;
	using OwnedSubMapView = typename Superclass::OwnedSubMapView;
	using AmountType = typename Superclass::AmountType;
	using ReactionDataRef = typename Superclass::ReactionDataRef;
	using ClusterData = typename Superclass::ClusterData;

	TrapReaction() = default;

	KOKKOS_INLINE_FUNCTION
	TrapReaction(ReactionDataRef reactionData, const ClusterData& clusterData,
		IndexType reactionId, IndexType cluster0, IndexType cluster1) :
		Superclass(reactionData, clusterData, reactionId),
		_reactant(cluster0),
		_trapped(cluster1)
	{
		this->initialize();
	}

	KOKKOS_INLINE_FUNCTION
	TrapReaction(ReactionDataRef reactionData, const ClusterData& clusterData,
		IndexType reactionId, const detail::ClusterSet& clusterSet) :
		TrapReaction(reactionData, clusterData, reactionId, clusterSet.cluster0,
			clusterSet.cluster1)
	{
	}

	static detail::CoefficientsView
	allocateCoefficientsView(IndexType)
	{
		return detail::CoefficientsView();
	}

	static detail::ConstantRateView
	allocateConstantRateView(IndexType, IndexType)
	{
		return detail::ConstantRateView();
	}

	//! Dummy method
	KOKKOS_FUNCTION
	double
	computeRate(IndexType gridIndex, double time = 0.0)
	{
		return 0.0;
	}

	KOKKOS_FUNCTION
	double
	computeTrappingRate(IndexType gridIndex, double time = 0.0);

	KOKKOS_FUNCTION
	double
	computeDeTrappingRate(IndexType gridIndex, double time = 0.0);

	KOKKOS_FUNCTION
	void
	setParameters(double den, double en, double str, double fre)
	{
		density = den;
		energy = en;
		strength = str;
		frequency = fre;
	}

	KOKKOS_FUNCTION
	IndexType
	getId()
	{
		return this->asDerived()->getId();
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
		// The reactant and trapped all connect together
		this->addConnectivity(_reactant, _reactant, connectivity);
		this->addConnectivity(_reactant, _trapped, connectivity);
		this->addConnectivity(_trapped, _reactant, connectivity);
		this->addConnectivity(_trapped, _trapped, connectivity);
	}

	KOKKOS_INLINE_FUNCTION
	void
	computeReducedConnectivity(const Connectivity& connectivity)
	{
		// The reactant connects with the reactant
		this->addConnectivity(_reactant, _reactant, connectivity);
		// The trapped connects with the trapped
		this->addConnectivity(_trapped, _trapped, connectivity);
	}

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
		IndexType gridIndex)
	{
		return 0.0;
	}

	KOKKOS_INLINE_FUNCTION
	void
	mapJacobianEntries(Connectivity connectivity)
	{
		_connEntries[0][0] = connectivity(_reactant, _reactant);
		_connEntries[0][1] = connectivity(_reactant, _trapped);
		_connEntries[1][0] = connectivity(_trapped, _reactant);
		_connEntries[1][1] = connectivity(_trapped, _trapped);
	}

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
	IndexType _trapped;
	static constexpr auto invalidIndex = Superclass::invalidIndex;

	util::Array<IndexType, 2, 2> _connEntries;

	// Trap parameters
	double density = 0.0;
	double energy = 0.0;
	double strength = 0.0;
	double frequency = 0.0;
};
} // namespace network
} // namespace core
} // namespace xolotl

#include <xolotl/core/network/detail/TrapReactionGenerator.h>
