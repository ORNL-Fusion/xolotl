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
SinkReactionGenerator<TBase>::SinkReactionGenerator(
	const NetworkType& network) :
	Superclass(network),
	_clusterSinkReactionCounts(
		"Sink Reaction Counts", Superclass::getNumberOfClusters())
{
}

template <typename TBase>
typename SinkReactionGenerator<TBase>::IndexType
SinkReactionGenerator<TBase>::getRowMapAndTotalReactionCount()
{
	_numPrecedingReactions = Superclass::getRowMapAndTotalReactionCount();
	_numSinkReactions = Kokkos::get_crs_row_map_from_counts(
		_sinkCrsRowMap, _clusterSinkReactionCounts);

	_sinkReactions =
		Kokkos::View<SinkReactionType*>("Sink Reactions", _numSinkReactions);

	return _numPrecedingReactions + _numSinkReactions;
}

template <typename TBase>
void
SinkReactionGenerator<TBase>::setupCrsClusterSetSubView()
{
	Superclass::setupCrsClusterSetSubView();
	_sinkCrsClusterSets = this->getClusterSetSubView(std::make_pair(
		_numPrecedingReactions, _numPrecedingReactions + _numSinkReactions));
}

template <typename TBase>
KOKKOS_INLINE_FUNCTION
void
SinkReactionGenerator<TBase>::addSinkReaction(
	Count, const ClusterSet& clusterSet) const
{
	if (!this->_clusterData.enableSink(0))
		return;

	Kokkos::atomic_increment(&_clusterSinkReactionCounts(clusterSet.cluster0));
}

template <typename TBase>
KOKKOS_INLINE_FUNCTION
void
SinkReactionGenerator<TBase>::addSinkReaction(
	Construct, const ClusterSet& clusterSet) const
{
	if (!this->_clusterData.enableSink(0))
		return;

	auto id = _sinkCrsRowMap(clusterSet.cluster0);
	for (; !Kokkos::atomic_compare_exchange_strong(
			 &_sinkCrsClusterSets(id).cluster0, NetworkType::invalidIndex(),
			 clusterSet.cluster0);
		 ++id) { }
	_sinkCrsClusterSets(id) = clusterSet;
}
} // namespace detail
} // namespace network
} // namespace core
} // namespace xolotl
