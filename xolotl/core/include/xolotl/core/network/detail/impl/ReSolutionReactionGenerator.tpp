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
ReSolutionReactionGenerator<TBase>::ReSolutionReactionGenerator(
	const NetworkType& network) :
	Superclass(network),
	_clusterReSoReactionCounts(
		"ReSolution Reaction Counts", Superclass::getNumberOfClusters())
{
}

template <typename TBase>
typename ReSolutionReactionGenerator<TBase>::IndexType
ReSolutionReactionGenerator<TBase>::getRowMapAndTotalReactionCount()
{
	_numPrecedingReactions = Superclass::getRowMapAndTotalReactionCount();
	_numReSoReactions = Kokkos::get_crs_row_map_from_counts(
		_reSoCrsRowMap, _clusterReSoReactionCounts);

	_reSoReactions = Kokkos::View<ReSolutionReactionType*>(
		"ReSolution Reactions", _numReSoReactions);

	return _numPrecedingReactions + _numReSoReactions;
}

template <typename TBase>
void
ReSolutionReactionGenerator<TBase>::setupCrsClusterSetSubView()
{
	Superclass::setupCrsClusterSetSubView();
	_reSoCrsClusterSets = this->getClusterSetSubView(std::make_pair(
		_numPrecedingReactions, _numPrecedingReactions + _numReSoReactions));
}

template <typename TBase>
KOKKOS_INLINE_FUNCTION
void
ReSolutionReactionGenerator<TBase>::addReSolutionReaction(
	Count, const ClusterSet& clusterSet) const
{
	if (!this->_clusterData.enableReSolution(0))
		return;

	Kokkos::atomic_increment(&_clusterReSoReactionCounts(clusterSet.cluster0));
}

template <typename TBase>
KOKKOS_INLINE_FUNCTION
void
ReSolutionReactionGenerator<TBase>::addReSolutionReaction(
	Construct, const ClusterSet& clusterSet) const
{
	if (!this->_clusterData.enableReSolution(0))
		return;

	auto id = _reSoCrsRowMap(clusterSet.cluster0);
	for (; !Kokkos::atomic_compare_exchange_strong(
			 &_reSoCrsClusterSets(id).cluster0, NetworkType::invalidIndex(),
			 clusterSet.cluster0);
		 ++id) { }
	_reSoCrsClusterSets(id) = clusterSet;
}
} // namespace detail
} // namespace network
} // namespace core
} // namespace xolotl
