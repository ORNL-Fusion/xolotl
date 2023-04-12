#pragma once

#include <type_traits>
#include <utility>

#include <xolotl/core/network/Reaction.h>
#include <xolotl/core/network/ReactionNetworkTraits.h>
#include <xolotl/core/network/detail/ReactionCollection.h>

namespace xolotl
{
namespace core
{
namespace network
{
namespace detail
{
/**
 * @brief General class filling collections for production
 * and dissociation reactions depending on the subpaving.
 *
 * @tparam TNetwork The network type
 * @tparam TDerived The derived class
 */
template <typename TNetwork, typename TDerived>
class ReactionGeneratorBase
{
public:
	using NetworkType = TNetwork;
	using NetworkTraits = ReactionNetworkTraits<NetworkType>;
	using ClusterData = typename NetworkType::ClusterData;
	using ClusterDataView = typename NetworkType::ClusterDataView;
	using Cluster = typename ClusterData::ClusterType;
	using ProductionReactionType =
		typename NetworkTraits::ProductionReactionType;
	using DissociationReactionType =
		typename NetworkTraits::DissociationReactionType;
	using Subpaving = typename NetworkType::Subpaving;
	using IndexType = typename NetworkType::IndexType;
	using IndexView = Kokkos::View<IndexType*>;
	using ClusterSetView = Kokkos::View<ClusterSet*>;
	using ClusterSetSubView =
		decltype(Kokkos::subview(std::declval<ClusterSetView>(),
			std::declval<std::pair<IndexType, IndexType>>()));
	using Connectivity = typename NetworkType::Connectivity;
	using ConnectivitiesView = typename NetworkType::ConnectivitiesView;

	struct Count
	{
	};

	struct Construct
	{
	};

	ReactionGeneratorBase(const TNetwork& network);

	ReactionCollection<NetworkType>
	generateReactions();

	KOKKOS_INLINE_FUNCTION
	const Subpaving&
	getSubpaving() const
	{
		return _subpaving;
	}

	KOKKOS_INLINE_FUNCTION
	Cluster
	getCluster(IndexType i) const
	{
		return _clusterData.getCluster(i);
	}

	IndexType
	getRowMapAndTotalReactionCount();

	ClusterSetSubView
	getClusterSetSubView(std::pair<IndexType, IndexType> indexRange)
	{
		return Kokkos::subview(_allClusterSets, indexRange);
	}

	void
	setupCrsClusterSetSubView();

	void
	setupCrs();

	KOKKOS_INLINE_FUNCTION
	IndexType
	getNumberOfClusters() const noexcept
	{
		return _clusterData.numClusters;
	}

	KOKKOS_INLINE_FUNCTION
	void
	addProductionReaction(Count, const ClusterSet& clusterSet) const;

	KOKKOS_INLINE_FUNCTION
	void
	addProductionReaction(Construct, const ClusterSet& clusterSet) const;

	KOKKOS_INLINE_FUNCTION
	void
	addDissociationReaction(Count, const ClusterSet& clusterSet) const;

	KOKKOS_INLINE_FUNCTION
	void
	addDissociationReaction(Construct, const ClusterSet& clusterSet) const;

	Kokkos::View<ProductionReactionType*>
	getProductionReactions() const
	{
		return _prodReactions;
	}

	Kokkos::View<DissociationReactionType*>
	getDissociationReactions() const
	{
		return _dissReactions;
	}

	void
	generateConnectivity(ReactionCollection<NetworkType>& reactionCollection);

	void
	setConstantConnectivities(ConnectivitiesView conns)
	{
		_constantConns = conns;
	}

	const ClusterConnectivity<>&
	getConnectivity() const
	{
		return _connectivity;
	}

protected:
	TDerived*
	asDerived()
	{
		return static_cast<TDerived*>(this);
	}

protected:
	Subpaving _subpaving;
	ClusterData _clusterData;
	ClusterDataView _clusterDataView;
	IndexType _numDOFs;
	bool _enableReducedJacobian;
	IndexView _clusterProdReactionCounts;
	IndexView _clusterDissReactionCounts;

	IndexType _numProdReactions;
	IndexType _numDissReactions;

	Kokkos::View<IndexType*> _prodCrsRowMap;
	Kokkos::View<IndexType*> _dissCrsRowMap;

	ClusterSetView _allClusterSets;
	ClusterSetSubView _prodCrsClusterSets;
	ClusterSetSubView _dissCrsClusterSets;

	Kokkos::View<ProductionReactionType*> _prodReactions;
	Kokkos::View<DissociationReactionType*> _dissReactions;

	ClusterConnectivity<> _connectivity;
	ConnectivitiesView _constantConns;
};

template <typename TNetwork, typename TReaction,
	typename TReactionGeneratorParent, typename = void>
struct WrapTypeSpecificReactionGenerator
{
};

template <typename TReactionGeneratorParent, typename TExtraReactionTypes>
struct ReactionGeneratorTypeBuilderImpl;

template <typename TReactionGeneratorParent>
struct ReactionGeneratorTypeBuilderImpl<TReactionGeneratorParent, std::tuple<>>
{
	using NetworkType = typename TReactionGeneratorParent::NetworkType;
	using Type = TReactionGeneratorParent;
};

template <typename TReactionGeneratorParent, typename... TExtraReactions>
struct ReactionGeneratorTypeBuilderImpl<TReactionGeneratorParent,
	std::tuple<TExtraReactions...>>
{
	using NetworkType = typename TReactionGeneratorParent::NetworkType;
	using ExtraReactions = std::tuple<TExtraReactions...>;
	using FrontReaction = std::tuple_element_t<0, ExtraReactions>;
	using Type =
		typename WrapTypeSpecificReactionGenerator<NetworkType, FrontReaction,
			typename ReactionGeneratorTypeBuilderImpl<TReactionGeneratorParent,
				TuplePopFront<ExtraReactions>>::Type>::Type;
};

template <typename TNetwork, typename TDerived>
struct ReactionGeneratorTypeBuilder
{
	using NetworkType = TNetwork;
	using ReactionTypes = ReactionTypeList<NetworkType>;
	using ReactionFirst = std::tuple_element_t<0, ReactionTypes>;
	using ReactionSecond = std::tuple_element_t<1, ReactionTypes>;

	static_assert(std::is_base_of_v<ProductionReaction<TNetwork, ReactionFirst>,
					  ReactionFirst>,
		"First reaction type must be a ProductionReaction");

	static_assert(
		std::is_base_of_v<DissociationReaction<TNetwork, ReactionSecond>,
			ReactionSecond>,
		"Second reaction type must be a DissociationReaction");

	using ExtraReactionTypes = TuplePopFront<TuplePopFront<ReactionTypes>>;
	using Type = typename ReactionGeneratorTypeBuilderImpl<
		ReactionGeneratorBase<TNetwork, TDerived>,
		TupleReverse<ExtraReactionTypes>>::Type;
};

template <typename TNetwork, typename TDerived>
using ReactionGenerator =
	typename ReactionGeneratorTypeBuilder<TNetwork, TDerived>::Type;
} // namespace detail
} // namespace network
} // namespace core
} // namespace xolotl
