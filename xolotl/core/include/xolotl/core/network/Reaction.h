#pragma once

#include <Kokkos_Core.hpp>

#include <plsm/Utility.h>
// #include <plsm/detail/KokkosExtension.h>

#include <xolotl/core/network/Cluster.h>
#include <xolotl/core/network/IReactionNetwork.h>
#include <xolotl/core/network/ReactionNetworkTraits.h>
#include <xolotl/core/network/SpeciesEnumSequence.h>
#include <xolotl/core/network/detail/ClusterSet.h>
#include <xolotl/core/network/detail/ReactionData.h>
#include <xolotl/util/Array.h>

namespace xolotl
{
namespace core
{
namespace network
{
/**
 * @brief General reaction class where
 * reactants become products with a given rate
 *
 * @tparam TNetwork The network type
 * @tparam TDerived The derived class type.
 */
template <typename TNetwork, typename TDerived>
class Reaction
{
	using Types = detail::ReactionNetworkTypes<TNetwork>;
	using Props = detail::ReactionNetworkProperties<TNetwork>;

protected:
	static constexpr auto invalidIndex = detail::invalidNetworkIndex;
	static constexpr auto nMomentIds = Props::numSpeciesNoI;
	static constexpr auto coeffsSingleExtent = Props::numSpeciesNoI + 1;

public:
	using NetworkType = TNetwork;
	using Species = typename Types::Species;
	using IndexType = typename Types::IndexType;
	using AmountType = typename Types::AmountType;
	using Region = typename Types::Region;
	using Composition = typename Types::Composition;
	using ConcentrationsView = IReactionNetwork::ConcentrationsView;
	using FluxesView = IReactionNetwork::FluxesView;
	using RatesView = IReactionNetwork::RatesView;
	using ConnectivitiesView = IReactionNetwork::ConnectivitiesView;
	using ConnectivitiesPairView = IReactionNetwork::ConnectivitiesPairView;
	using BelongingView = IReactionNetwork::BelongingView;
	using OwnedSubMapView = IReactionNetwork::OwnedSubMapView;
	using Connectivity = typename IReactionNetwork::Connectivity;
	using ReactionDataRef = typename Types::ReactionDataRef;
	using ClusterData = typename Types::ClusterData;
	using ReflectedRegion =
		plsm::Region<plsm::DifferenceType<typename Region::ScalarType>,
			Props::numSpeciesNoI>;

	Reaction() = default;

	KOKKOS_INLINE_FUNCTION
	Reaction(ReactionDataRef reactionData, const ClusterData& clusterData,
		IndexType reactionId);

	KOKKOS_INLINE_FUNCTION
	void
	updateData(ReactionDataRef reactionData, const ClusterData& clusterData);

	KOKKOS_INLINE_FUNCTION
	void
	getRateEntries(ReactionDataRef reactionData);

	KOKKOS_INLINE_FUNCTION
	void
	updateRates(double time = 0.0)
	{
		for (IndexType i = 0; i < _rate.extent(0); ++i) {
			_rate(i) = asDerived()->computeRate(i, time);
		}
	}

	/**
	 * @brief Computes the contribution to the connectivity
	 * (which cluster interacts with which one).
	 */
	KOKKOS_INLINE_FUNCTION
	void
	contributeConnectivity(const Connectivity& connectivity)
	{
		asDerived()->computeConnectivity(connectivity);
	}

	/**
	 * @brief Computes the contribution to the connectivity
	 * when the reduced matrix method is used (only the
	 * diagonal)
	 */
	KOKKOS_INLINE_FUNCTION
	void
	contributeReducedConnectivity(const Connectivity& connectivity)
	{
		asDerived()->computeReducedConnectivity(connectivity);
	}

	KOKKOS_INLINE_FUNCTION
	void
	contributeFlux(ConcentrationsView concentrations, FluxesView fluxes,
		IndexType gridIndex)
	{
		asDerived()->computeFlux(concentrations, fluxes, gridIndex);
	}

	/**
	 * @brief Computes the contribution to the Jacobian.
	 */
	KOKKOS_INLINE_FUNCTION
	void
	contributePartialDerivatives(ConcentrationsView concentrations,
		Kokkos::View<double*> values, IndexType gridIndex)
	{
		asDerived()->computePartialDerivatives(
			concentrations, values, gridIndex);
	}

	/**
	 * @brief Computes the contribution to the Jacobian
	 * when the reduced matrix method is used (only the
	 * diagonal)
	 */
	KOKKOS_INLINE_FUNCTION
	void
	contributeReducedPartialDerivatives(ConcentrationsView concentrations,
		Kokkos::View<double*> values, IndexType gridIndex)
	{
		asDerived()->computeReducedPartialDerivatives(
			concentrations, values, gridIndex);
	}

	KOKKOS_INLINE_FUNCTION
	void
	contributeConstantRates(ConcentrationsView concentrations, RatesView rates,
		BelongingView isInSub, IndexType subId, IndexType gridIndex)
	{
		asDerived()->computeConstantRates(
			concentrations, rates, isInSub, subId, gridIndex);
	}

	KOKKOS_INLINE_FUNCTION
	void
	contributeConstantConnectivities(ConnectivitiesView conns,
		BelongingView isInSub, OwnedSubMapView backMap)
	{
		asDerived()->getConstantConnectivities(conns, isInSub, backMap);
	}

	/**
	 * \see IReactionNetwork.h
	 */
	KOKKOS_INLINE_FUNCTION
	double
	contributeLeftSideRate(ConcentrationsView concentrations,
		IndexType clusterId, IndexType gridIndex)
	{
		return asDerived()->computeLeftSideRate(
			concentrations, clusterId, gridIndex);
	}

	KOKKOS_INLINE_FUNCTION
	void
	defineJacobianEntries(Connectivity connectivity)
	{
		asDerived()->mapJacobianEntries(connectivity);
	}

	KOKKOS_INLINE_FUNCTION
	void
	defineRateEntries(ConnectivitiesPairView connectivityRow,
		ConnectivitiesPairView connectivityEntries,
		BelongingView isInSub = BelongingView(),
		OwnedSubMapView backMap = OwnedSubMapView(), IndexType subId = 0)
	{
		asDerived()->mapRateEntries(
			connectivityRow, connectivityEntries, isInSub, backMap, subId);
	}

protected:
	KOKKOS_INLINE_FUNCTION
	TDerived*
	asDerived()
	{
		return static_cast<TDerived*>(this);
	}

	KOKKOS_INLINE_FUNCTION
	void
	initialize()
	{
		asDerived()->computeCoefficients();
		updateRates();
	}

	/**
	 * @brief Computes the volume by which the reactants and products
	 * overlap, making the reaction viable.
	 */
	KOKKOS_INLINE_FUNCTION
	double
	computeOverlap(const ReflectedRegion& cl1RR, const ReflectedRegion& cl2RR,
		const ReflectedRegion& pr1RR, const ReflectedRegion& pr2RR);

	KOKKOS_INLINE_FUNCTION
	void
	copyMomentIds(
		IndexType clusterId, util::Array<IndexType, nMomentIds>& momentIds)
	{
		if (clusterId == invalidIndex) {
			for (IndexType i = 0; i < nMomentIds; ++i) {
				momentIds[i] = invalidIndex;
			}
			return;
		}

		const auto& mIds = _clusterData->getCluster(clusterId).getMomentIds();
		for (IndexType i = 0; i < nMomentIds; ++i) {
			momentIds[i] = mIds[i];
		}
	}

	KOKKOS_INLINE_FUNCTION
	void
	addConnectivity(
		IndexType rowId, IndexType columnId, const Connectivity& connectivity)
	{
		connectivity.add(rowId, columnId);
	}

	KOKKOS_INLINE_FUNCTION
	IndexType
	getPosition(IndexType rowId, IndexType columnId,
		const ConnectivitiesPairView connectivityRow,
		const ConnectivitiesPairView connectivityEntries) const
	{
		for (auto pos = connectivityRow(rowId);
			 pos < connectivityRow(rowId + 1); ++pos) {
			if (connectivityEntries(pos) == columnId) {
				return pos;
			}
		}
		return invalidIndex;
	}

protected:
	const ClusterData* _clusterData;

	IndexType _reactionId{invalidIndex};
	double _deltaG0;

	//! Reaction rate (k)
	using RateSubView = decltype(std::declval<ReactionDataRef>().getRates(0));
	RateSubView _rate;

	//! Reaction widths
	using WidthSubView = decltype(std::declval<ReactionDataRef>().getWidths(0));
	WidthSubView _widths;

	//! Flux coefficients
	using CoefsSubView =
		decltype(std::declval<ReactionDataRef>().getCoefficients(0));
	CoefsSubView _coefs;

	//! Constant Rates (only actually used for constant reactions
	using ConstantRateSubView =
		decltype(std::declval<ReactionDataRef>().getConstantRates(0));
	ConstantRateSubView _constantRates;

	using RateEntriesSubView =
		decltype(std::declval<ReactionDataRef>().getRateEntries(0));
	RateEntriesSubView _rateEntries;
};

/**
 * @brief Class implementing production reaction where
 * R_1 + R_2 -> sum_n P_n with a given rate, and n = 0, 1, 2.
 *
 * @tparam TNetwork The network type
 * @tparam TDerived The derived class type.
 */
template <typename TNetwork, typename TDerived>
class ProductionReaction : public Reaction<TNetwork, TDerived>
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
	using Composition = typename Superclass::Composition;
	using Region = typename Superclass::Region;
	using AmountType = typename Superclass::AmountType;
	using ReactionDataRef = typename Superclass::ReactionDataRef;
	using ClusterData = typename Superclass::ClusterData;
	using ReflectedRegion = typename Superclass::ReflectedRegion;

	ProductionReaction() = default;

	KOKKOS_INLINE_FUNCTION
	ProductionReaction(ReactionDataRef reactionData,
		const ClusterData& clusterData, IndexType reactionId,
		IndexType cluster0, IndexType cluster1,
		IndexType cluster2 = invalidIndex, IndexType cluster3 = invalidIndex);

	KOKKOS_INLINE_FUNCTION
	ProductionReaction(ReactionDataRef reactionData,
		const ClusterData& clusterData, IndexType reactionId,
		const detail::ClusterSet& clusterSet);

	static detail::CoefficientsView
	allocateCoefficientsView(IndexType size)
	{
		return detail::CoefficientsView("Production Coefficients", size,
			Superclass::coeffsSingleExtent, Superclass::coeffsSingleExtent, 4,
			Superclass::coeffsSingleExtent);
	}

	static detail::ConstantRateView
	allocateConstantRateView(IndexType, IndexType)
	{
		return detail::ConstantRateView();
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
		BelongingView isInSub, IndexType subId, IndexType gridIndex);

	KOKKOS_INLINE_FUNCTION
	void
	getConstantConnectivities(ConnectivitiesView conns, BelongingView isInSub,
		OwnedSubMapView backMap);

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
		OwnedSubMapView backMap, IndexType subId);

protected:
	static constexpr auto invalidIndex = Superclass::invalidIndex;
	util::Array<IndexType, 2> _reactants{invalidIndex, invalidIndex};
	util::Array<IndexType, 2> _products{invalidIndex, invalidIndex};
	util::Array<double, 2> _reactantVolumes{0.0, 0.0};
	util::Array<double, 2> _productVolumes{0.0, 0.0};

	static constexpr auto nMomentIds = Superclass::nMomentIds;
	util::Array<IndexType, 2, nMomentIds> _reactantMomentIds;
	util::Array<IndexType, 2, nMomentIds> _productMomentIds;

	util::Array<IndexType, 4, 1 + nMomentIds, 2, 1 + nMomentIds> _connEntries;
};

/**
 * @brief Class implementing dissociation reaction where
 * R -> P_1 + P_2 with a given rate.
 *
 * @tparam TNetwork The network type
 * @tparam TDerived The derived class type.
 */
template <typename TNetwork, typename TDerived>
class DissociationReaction : public Reaction<TNetwork, TDerived>
{
	friend class Reaction<TNetwork, TDerived>;

public:
	using NetworkType = TNetwork;
	using Superclass = Reaction<TNetwork, TDerived>;
	using IndexType = typename Superclass::IndexType;
	using Region = typename Superclass::Region;
	using Composition = typename Superclass::Composition;
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
	using ReflectedRegion = typename Superclass::ReflectedRegion;

	DissociationReaction() = default;

	KOKKOS_INLINE_FUNCTION
	DissociationReaction(ReactionDataRef reactionData,
		const ClusterData& clusterData, IndexType reactionId,
		IndexType cluster0, IndexType cluster1, IndexType cluster2);

	KOKKOS_INLINE_FUNCTION
	DissociationReaction(ReactionDataRef reactionData,
		const ClusterData& clusterData, IndexType reactionId,
		const detail::ClusterSet& clusterSet);

	static detail::CoefficientsView
	allocateCoefficientsView(IndexType size)
	{
		return detail::CoefficientsView("Dissociation Coefficients", size,
			Superclass::coeffsSingleExtent, 1, 3,
			Superclass::coeffsSingleExtent);
	}

	static detail::ConstantRateView
	allocateConstantRateView(IndexType, IndexType)
	{
		return detail::ConstantRateView();
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
		BelongingView isInSub, IndexType subId, IndexType gridIndex);

	KOKKOS_INLINE_FUNCTION
	void
	getConstantConnectivities(ConnectivitiesView conns, BelongingView isInSub,
		OwnedSubMapView backMap);

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
		OwnedSubMapView backMap, IndexType subId);

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
