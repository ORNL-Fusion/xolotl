#pragma once

namespace xolotl
{
namespace core
{
namespace network
{
namespace detail
{
template <typename TBase>
TrapMutationReactionGenerator<TBase>::TrapMutationReactionGenerator(
	const NetworkType& network) :
	Superclass(network),
	_clusterTMReactionCounts(
		"Trap Mutation Reaction Counts", Superclass::getNumberOfClusters())
{
}

template <typename TBase>
typename TrapMutationReactionGenerator<TBase>::IndexType
TrapMutationReactionGenerator<TBase>::getRowMapAndTotalReactionCount()
{
	_numPrecedingReactions = Superclass::getRowMapAndTotalReactionCount();
	_numTMReactions = Kokkos::get_crs_row_map_from_counts(
		_tmCrsRowMap, _clusterTMReactionCounts);

	_tmReactions = Kokkos::View<TrapMutationReactionType*>(
		"Trap Mutation Reactions", _numTMReactions);

	return _numPrecedingReactions + _numTMReactions;
}

template <typename TBase>
void
TrapMutationReactionGenerator<TBase>::setupCrsClusterSetSubView()
{
	Superclass::setupCrsClusterSetSubView();
	_tmCrsClusterSets = this->getClusterSetSubView(std::make_pair(
		_numPrecedingReactions, _numPrecedingReactions + _numTMReactions));
}

template <typename TBase>
KOKKOS_INLINE_FUNCTION
void
TrapMutationReactionGenerator<TBase>::addTrapMutationReaction(
	Count, const ClusterSet& clusterSet) const
{
	Kokkos::atomic_increment(&_clusterTMReactionCounts(clusterSet.cluster1));
}

template <typename TBase>
KOKKOS_INLINE_FUNCTION
void
TrapMutationReactionGenerator<TBase>::addTrapMutationReaction(
	Construct, const ClusterSet& clusterSet) const
{
	auto id = _tmCrsRowMap(clusterSet.cluster1);
	for (; !Kokkos::atomic_compare_exchange_strong(
			 &_tmCrsClusterSets(id).cluster1, NetworkType::invalidIndex(),
			 clusterSet.cluster1);
		 ++id) { }
	_tmCrsClusterSets(id) = clusterSet;
}
} // namespace detail
} // namespace network
} // namespace core
} // namespace xolotl
