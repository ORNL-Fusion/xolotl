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
TrapReactionGenerator<TBase>::TrapReactionGenerator(
	const NetworkType& network) :
	Superclass(network),
	_clusterTrapReactionCounts(
		"Trap Reaction Counts", Superclass::getNumberOfClusters())
{
}

template <typename TBase>
typename TrapReactionGenerator<TBase>::IndexType
TrapReactionGenerator<TBase>::getRowMapAndTotalReactionCount()
{
	_numPrecedingReactions = Superclass::getRowMapAndTotalReactionCount();
	_numTrapReactions = Kokkos::get_crs_row_map_from_counts(
		_trapCrsRowMap, _clusterTrapReactionCounts);

	_trapReactions =
		Kokkos::View<TrapReactionType*>("Trap Reactions", _numTrapReactions);

	return _numPrecedingReactions + _numTrapReactions;
}

template <typename TBase>
void
TrapReactionGenerator<TBase>::setupCrsClusterSetSubView()
{
	Superclass::setupCrsClusterSetSubView();
	_trapCrsClusterSets = this->getClusterSetSubView(std::make_pair(
		_numPrecedingReactions, _numPrecedingReactions + _numTrapReactions));
}

template <typename TBase>
KOKKOS_INLINE_FUNCTION
void
TrapReactionGenerator<TBase>::addTrapReaction(
	Count, const ClusterSet& clusterSet) const
{
	if (!this->_clusterData.enableTrap())
		return;

	Kokkos::atomic_inc(&_clusterTrapReactionCounts(clusterSet.cluster1));
}

template <typename TBase>
KOKKOS_INLINE_FUNCTION
void
TrapReactionGenerator<TBase>::addTrapReaction(
	Construct, const ClusterSet& clusterSet) const
{
	if (!this->_clusterData.enableTrap())
		return;

	auto id = _trapCrsRowMap(clusterSet.cluster0);
	for (; !util::atomicCompareExchangeStrong(&_trapCrsClusterSets(id).cluster1,
			 NetworkType::invalidIndex(), clusterSet.cluster0);
		++id) { }
	_trapCrsClusterSets(id) = clusterSet;
}
} // namespace detail
} // namespace network
} // namespace core
} // namespace xolotl
