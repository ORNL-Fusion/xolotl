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
PartialReSolutionReactionGenerator<TBase>::PartialReSolutionReactionGenerator(
	const NetworkType& network) :
	Superclass(network),
	_clusterPartialReSoReactionCounts(
		"PartialReSolution Reaction Counts", Superclass::getNumberOfClusters())
{
}

template <typename TBase>
typename PartialReSolutionReactionGenerator<TBase>::IndexType
PartialReSolutionReactionGenerator<TBase>::getRowMapAndTotalReactionCount()
{
	_numPrecedingReactions = Superclass::getRowMapAndTotalReactionCount();
	_numPartialReSoReactions = Kokkos::get_crs_row_map_from_counts(
		_reSoCrsRowMap, _clusterPartialReSoReactionCounts);

	_reSoReactions = Kokkos::View<PartialReSolutionReactionType*>(
		"PartialReSolution Reactions", _numPartialReSoReactions);

	return _numPrecedingReactions + _numPartialReSoReactions;
}

template <typename TBase>
void
PartialReSolutionReactionGenerator<TBase>::setupCrsClusterSetSubView()
{
	Superclass::setupCrsClusterSetSubView();
	_reSoCrsClusterSets =
		this->getClusterSetSubView(std::make_pair(_numPrecedingReactions,
			_numPrecedingReactions + _numPartialReSoReactions));
}

template <typename TBase>
KOKKOS_INLINE_FUNCTION
void
PartialReSolutionReactionGenerator<TBase>::addPartialReSolutionReaction(
	Count, const ClusterSet& clusterSet) const
{
	if (!this->_clusterData.enablePartialReSolution())
		return;

	Kokkos::atomic_increment(
		&_clusterPartialReSoReactionCounts(clusterSet.cluster0));
}

template <typename TBase>
KOKKOS_INLINE_FUNCTION
void
PartialReSolutionReactionGenerator<TBase>::addPartialReSolutionReaction(
	Construct, const ClusterSet& clusterSet) const
{
	if (!this->_clusterData.enablePartialReSolution())
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
