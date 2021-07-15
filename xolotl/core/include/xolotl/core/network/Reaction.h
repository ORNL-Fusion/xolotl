#pragma once

#include <Kokkos_Core.hpp>

#include <plsm/Utility.h>
#include <plsm/detail/KokkosExtension.h>

#include <xolotl/core/network/Cluster.h>
#include <xolotl/core/network/IReactionNetwork.h>
#include <xolotl/core/network/ReactionNetworkTraits.h>
#include <xolotl/core/network/SpeciesEnumSequence.h>
#include <xolotl/core/network/detail/ClusterSet.h>
#include <xolotl/core/network/detail/ReactionData.h>

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
	using ClusterDataRef = typename Types::ClusterDataRef;
	using ConcentrationsView = IReactionNetwork::ConcentrationsView;
	using FluxesView = IReactionNetwork::FluxesView;
	using Connectivity = typename IReactionNetwork::Connectivity;
	using ReactionDataRef = typename Types::ReactionDataRef;
	using ReflectedRegion =
		plsm::Region<plsm::DifferenceType<typename Region::ScalarType>,
			Props::numSpeciesNoI>;

	Reaction() = default;

	KOKKOS_INLINE_FUNCTION
	Reaction(ReactionDataRef reactionData, ClusterDataRef clusterData,
		IndexType reactionId);

	KOKKOS_INLINE_FUNCTION
	void
	updateData(ReactionDataRef reactionData, ClusterDataRef clusterData);

	KOKKOS_INLINE_FUNCTION
	void
	updateRates()
	{
		for (IndexType i = 0; i < _rate.extent(0); ++i) {
			_rate(i) = asDerived()->computeRate(i);
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
		Kokkos::View<double*> values, Connectivity connectivity,
		IndexType gridIndex)
	{
		asDerived()->computePartialDerivatives(
			concentrations, values, connectivity, gridIndex);
	}

	/**
	 * @brief Computes the contribution to the Jacobian
	 * when the reduced matrix method is used (only the
	 * diagonal)
	 */
	KOKKOS_INLINE_FUNCTION
	void
	contributeReducedPartialDerivatives(ConcentrationsView concentrations,
		Kokkos::View<double*> values, Connectivity connectivity,
		IndexType gridIndex)
	{
		asDerived()->computeReducedPartialDerivatives(
			concentrations, values, connectivity, gridIndex);
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
	AmountType
	computeOverlap(const ReflectedRegion& cl1RR, const ReflectedRegion& cl2RR,
		const ReflectedRegion& pr1RR, const ReflectedRegion& pr2RR);

	KOKKOS_INLINE_FUNCTION
	void
	copyMomentIds(
		IndexType clusterId, Kokkos::Array<IndexType, nMomentIds>& momentIds)
	{
		if (clusterId == invalidIndex) {
			for (IndexType i = 0; i < nMomentIds; ++i) {
				momentIds[i] = invalidIndex;
			}
			return;
		}

		const auto& mIds = _clusterData.getCluster(clusterId).getMomentIds();
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

protected:
	ClusterDataRef _clusterData;

	IndexType _reactionId{invalidIndex};

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
	using ClusterDataRef = typename Superclass::ClusterDataRef;
	using IndexType = typename Superclass::IndexType;
	using Connectivity = typename Superclass::Connectivity;
	using ConcentrationsView = typename Superclass::ConcentrationsView;
	using FluxesView = typename Superclass::FluxesView;
	using Composition = typename Superclass::Composition;
	using Region = typename Superclass::Region;
	using AmountType = typename Superclass::AmountType;
	using ReactionDataRef = typename Superclass::ReactionDataRef;
	using ReflectedRegion = typename Superclass::ReflectedRegion;

	ProductionReaction() = default;

	KOKKOS_INLINE_FUNCTION
	ProductionReaction(ReactionDataRef reactionData, ClusterDataRef clusterData,
		IndexType reactionId, IndexType cluster0, IndexType cluster1,
		IndexType cluster2 = invalidIndex, IndexType cluster3 = invalidIndex);

	KOKKOS_INLINE_FUNCTION
	ProductionReaction(ReactionDataRef reactionData, ClusterDataRef clusterData,
		IndexType reactionId, const detail::ClusterSet& clusterSet);

	static detail::CoefficientsView
	allocateCoefficientsView(IndexType size)
	{
		return detail::CoefficientsView("Production Coefficients", size,
			Superclass::coeffsSingleExtent, Superclass::coeffsSingleExtent, 4,
			Superclass::coeffsSingleExtent);
	}

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
	static constexpr auto invalidIndex = Superclass::invalidIndex;
	Kokkos::Array<IndexType, 2> _reactants{invalidIndex, invalidIndex};
	Kokkos::Array<IndexType, 2> _products{invalidIndex, invalidIndex};
	Kokkos::Array<AmountType, 2> _reactantVolumes{0, 0};
	Kokkos::Array<AmountType, 2> _productVolumes{0, 0};

	static constexpr auto nMomentIds = Superclass::nMomentIds;
	Kokkos::Array<Kokkos::Array<IndexType, nMomentIds>, 2> _reactantMomentIds;
	Kokkos::Array<Kokkos::Array<IndexType, nMomentIds>, 2> _productMomentIds;
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
	using ClusterDataRef = typename Superclass::ClusterDataRef;
	using IndexType = typename Superclass::IndexType;
	using Region = typename Superclass::Region;
	using Composition = typename Superclass::Composition;
	using Connectivity = typename Superclass::Connectivity;
	using ConcentrationsView = typename Superclass::ConcentrationsView;
	using FluxesView = typename Superclass::FluxesView;
	using AmountType = typename Superclass::AmountType;
	using ReactionDataRef = typename Superclass::ReactionDataRef;
	using ReflectedRegion = typename Superclass::ReflectedRegion;

	DissociationReaction() = default;

	KOKKOS_INLINE_FUNCTION
	DissociationReaction(ReactionDataRef reactionData,
		ClusterDataRef clusterData, IndexType reactionId, IndexType cluster0,
		IndexType cluster1, IndexType cluster2);

	KOKKOS_INLINE_FUNCTION
	DissociationReaction(ReactionDataRef reactionData,
		ClusterDataRef clusterData, IndexType reactionId,
		const detail::ClusterSet& clusterSet);

	static detail::CoefficientsView
	allocateCoefficientsView(IndexType size)
	{
		return detail::CoefficientsView("Dissociation Coefficients", size,
			Superclass::coeffsSingleExtent, 1, 3,
			Superclass::coeffsSingleExtent);
	}

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
	AmountType _reactantVolume;
	static constexpr auto invalidIndex = Superclass::invalidIndex;
	Kokkos::Array<IndexType, 2> _products{invalidIndex, invalidIndex};
	Kokkos::Array<AmountType, 2> _productVolumes{0, 0};

	static constexpr auto nMomentIds = Superclass::nMomentIds;
	Kokkos::Array<IndexType, nMomentIds> _reactantMomentIds;
	Kokkos::Array<Kokkos::Array<IndexType, nMomentIds>, 2> _productMomentIds;
};
} // namespace network
} // namespace core
} // namespace xolotl
