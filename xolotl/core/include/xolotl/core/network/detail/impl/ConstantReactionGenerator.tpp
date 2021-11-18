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
ConstantReactionGenerator<TBase>::ConstantReactionGenerator(
	const NetworkType& network) :
	Superclass(network),
	_clusterConstantReactionCounts(
		"Constant Reaction Counts", Superclass::getNumberOfClusters())
{
}

template <typename TBase>
typename ConstantReactionGenerator<TBase>::IndexType
ConstantReactionGenerator<TBase>::getRowMapAndTotalReactionCount()
{
	_numPrecedingReactions = Superclass::getRowMapAndTotalReactionCount();
	_numConstantReactions = Kokkos::get_crs_row_map_from_counts(
		_constantCrsRowMap, _clusterConstantReactionCounts);

	_constantReactions = Kokkos::View<ConstantReactionType*>(
		"Constant Reactions", _numConstantReactions);

	return _numPrecedingReactions + _numConstantReactions;
}

template <typename TBase>
void
ConstantReactionGenerator<TBase>::setupCrsClusterSetSubView()
{
	Superclass::setupCrsClusterSetSubView();
	_constantCrsClusterSets =
		this->getClusterSetSubView(std::make_pair(_numPrecedingReactions,
			_numPrecedingReactions + _numConstantReactions));
}

template <typename TBase>
KOKKOS_INLINE_FUNCTION
void
ConstantReactionGenerator<TBase>::addConstantReaction(
	Count, const ClusterSet& clusterSet) const
{
	if (!this->_clusterData.enableConstantReaction())
		return;

	Kokkos::atomic_increment(
		&_clusterConstantReactionCounts(clusterSet.cluster0));
}

template <typename TBase>
KOKKOS_INLINE_FUNCTION
void
ConstantReactionGenerator<TBase>::addConstantReaction(
	Construct, const ClusterSet& clusterSet) const
{
	if (!this->_clusterData.enableConstantReaction())
		return;

	auto id = _constantCrsRowMap(clusterSet.cluster0);
	for (; !Kokkos::atomic_compare_exchange_strong(
			 &_constantCrsClusterSets(id).cluster0, NetworkType::invalidIndex(),
			 clusterSet.cluster0);
		 ++id) { }
	_constantCrsClusterSets(id) = clusterSet;
}
} // namespace detail
} // namespace network
} // namespace core
} // namespace xolotl
