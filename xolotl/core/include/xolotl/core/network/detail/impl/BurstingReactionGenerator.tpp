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
BurstingReactionGenerator<TBase>::BurstingReactionGenerator(
	const NetworkType& network) :
	Superclass(network),
	_clusterBurstingReactionCounts(
		"Bursting Reaction Counts", Superclass::getNumberOfClusters())
{
}

template <typename TBase>
typename BurstingReactionGenerator<TBase>::IndexType
BurstingReactionGenerator<TBase>::getRowMapAndTotalReactionCount()
{
	_numPrecedingReactions = Superclass::getRowMapAndTotalReactionCount();
	_numBurstingReactions = Kokkos::get_crs_row_map_from_counts(
		_burstingCrsRowMap, _clusterBurstingReactionCounts);

	_burstingReactions = Kokkos::View<BurstingReactionType*>(
		"Bursting Reactions", _numBurstingReactions);

	return _numPrecedingReactions + _numBurstingReactions;
}

template <typename TBase>
void
BurstingReactionGenerator<TBase>::setupCrsClusterSetSubView()
{
	Superclass::setupCrsClusterSetSubView();
	_burstingCrsClusterSets =
		this->getClusterSetSubView(std::make_pair(_numPrecedingReactions,
			_numPrecedingReactions + _numBurstingReactions));
}

template <typename TBase>
KOKKOS_INLINE_FUNCTION
void
BurstingReactionGenerator<TBase>::addBurstingReaction(
	Count, const ClusterSet& clusterSet) const
{
	if (!this->_clusterData.enableBurst())
		return;

	Kokkos::atomic_increment(
		&_clusterBurstingReactionCounts(clusterSet.cluster0));
}

template <typename TBase>
KOKKOS_INLINE_FUNCTION
void
BurstingReactionGenerator<TBase>::addBurstingReaction(
	Construct, const ClusterSet& clusterSet) const
{
	if (!this->_clusterData.enableBurst())
		return;

	auto id = _burstingCrsRowMap(clusterSet.cluster0);
	for (; !Kokkos::atomic_compare_exchange_strong(
			 &_burstingCrsClusterSets(id).cluster0, NetworkType::invalidIndex(),
			 clusterSet.cluster0);
		 ++id) { }
	_burstingCrsClusterSets(id) = clusterSet;
}
} // namespace detail
} // namespace network
} // namespace core
} // namespace xolotl
