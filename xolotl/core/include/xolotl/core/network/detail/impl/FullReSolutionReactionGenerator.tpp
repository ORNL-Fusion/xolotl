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
FullReSolutionReactionGenerator<TBase>::FullReSolutionReactionGenerator(
	const NetworkType& network) :
	Superclass(network),
	_clusterFullReSoReactionCounts(
		"FullReSolution Reaction Counts", Superclass::getNumberOfClusters())
{
}

template <typename TBase>
typename FullReSolutionReactionGenerator<TBase>::IndexType
FullReSolutionReactionGenerator<TBase>::getRowMapAndTotalReactionCount()
{
	_numPrecedingReactions = Superclass::getRowMapAndTotalReactionCount();
	_numFullReSoReactions = Kokkos::get_crs_row_map_from_counts(
		_reSoCrsRowMap, _clusterFullReSoReactionCounts);

	_reSoReactions = Kokkos::View<FullReSolutionReactionType*>(
		"FullReSolution Reactions", _numFullReSoReactions);

	return _numPrecedingReactions + _numFullReSoReactions;
}

template <typename TBase>
void
FullReSolutionReactionGenerator<TBase>::setupCrsClusterSetSubView()
{
	Superclass::setupCrsClusterSetSubView();
	_reSoCrsClusterSets =
		this->getClusterSetSubView(std::make_pair(_numPrecedingReactions,
			_numPrecedingReactions + _numFullReSoReactions));
}

template <typename TBase>
KOKKOS_INLINE_FUNCTION
void
FullReSolutionReactionGenerator<TBase>::addFullReSolutionReaction(
	Count, const ClusterSet& clusterSet) const
{
	if (!this->_clusterData.enableFullReSolution())
		return;

	Kokkos::atomic_increment(
		&_clusterFullReSoReactionCounts(clusterSet.cluster0));
}

template <typename TBase>
KOKKOS_INLINE_FUNCTION
void
FullReSolutionReactionGenerator<TBase>::addFullReSolutionReaction(
	Construct, const ClusterSet& clusterSet) const
{
	if (!this->_clusterData.enableFullReSolution())
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
