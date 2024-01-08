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
TransformReactionGenerator<TBase>::TransformReactionGenerator(
	const NetworkType& network) :
	Superclass(network),
	_clusterTransformReactionCounts(
		"Transform Reaction Counts", Superclass::getNumberOfClusters())
{
}

template <typename TBase>
typename TransformReactionGenerator<TBase>::IndexType
TransformReactionGenerator<TBase>::getRowMapAndTotalReactionCount()
{
	_numPrecedingReactions = Superclass::getRowMapAndTotalReactionCount();
	_numTransformReactions = Kokkos::get_crs_row_map_from_counts(
		_transformCrsRowMap, _clusterTransformReactionCounts);

	_transformReactions = Kokkos::View<TransformReactionType*>(
		"Transform Reactions", _numTransformReactions);

	return _numPrecedingReactions + _numTransformReactions;
}

template <typename TBase>
void
TransformReactionGenerator<TBase>::setupCrsClusterSetSubView()
{
	Superclass::setupCrsClusterSetSubView();
	_transformCrsClusterSets =
		this->getClusterSetSubView(std::make_pair(_numPrecedingReactions,
			_numPrecedingReactions + _numTransformReactions));
}

template <typename TBase>
KOKKOS_INLINE_FUNCTION
void
TransformReactionGenerator<TBase>::addTransformReaction(
	Count, const ClusterSet& clusterSet) const
{
	if (!this->_clusterData.enableStdReaction())
		return;

	Kokkos::atomic_increment(
		&_clusterTransformReactionCounts(clusterSet.cluster0));
}

template <typename TBase>
KOKKOS_INLINE_FUNCTION
void
TransformReactionGenerator<TBase>::addTransformReaction(
	Construct, const ClusterSet& clusterSet) const
{
	if (!this->_clusterData.enableStdReaction())
		return;

	auto id = _transformCrsRowMap(clusterSet.cluster0);
	for (; !Kokkos::atomic_compare_exchange_strong(
			 &_transformCrsClusterSets(id).cluster0,
			 NetworkType::invalidIndex(), clusterSet.cluster0);
		 ++id) { }
	_transformCrsClusterSets(id) = clusterSet;
}
} // namespace detail
} // namespace network
} // namespace core
} // namespace xolotl
