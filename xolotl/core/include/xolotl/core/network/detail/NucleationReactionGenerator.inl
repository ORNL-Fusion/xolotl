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
NucleationReactionGenerator<TBase>::NucleationReactionGenerator(
        const NetworkType& network)
    :
    Superclass(network),
    _clusterNucleationReactionCounts("Nucleation Reaction Counts",
        Superclass::getNumberOfClusters())
{
}

template <typename TBase>
typename NucleationReactionGenerator<TBase>::IndexType
NucleationReactionGenerator<TBase>::getRowMapAndTotalReactionCount()
{
    _numPrecedingReactions = Superclass::getRowMapAndTotalReactionCount();
    _numNucleationReactions = Kokkos::get_crs_row_map_from_counts(_nucleationCrsRowMap,
        _clusterNucleationReactionCounts);

    _nucleationReactions = Kokkos::View<NucleationReactionType*>("Nucleation Reactions",
        _numNucleationReactions);

    return _numPrecedingReactions + _numNucleationReactions;
}

template <typename TBase>
void
NucleationReactionGenerator<TBase>::setupCrsClusterSetSubView()
{
    Superclass::setupCrsClusterSetSubView();
    _nucleationCrsClusterSets = this->getClusterSetSubView(
        std::make_pair(_numPrecedingReactions,
            _numPrecedingReactions + _numNucleationReactions));
}

template <typename TBase>
KOKKOS_INLINE_FUNCTION
void
NucleationReactionGenerator<TBase>::addNucleationReaction(Count,
    const ClusterSet& clusterSet) const
{
    Kokkos::atomic_increment(
        &_clusterNucleationReactionCounts(clusterSet.cluster0));
}

template <typename TBase>
KOKKOS_INLINE_FUNCTION
void
NucleationReactionGenerator<TBase>::addNucleationReaction(Construct,
    const ClusterSet& clusterSet) const
{
    auto id = _nucleationCrsRowMap(clusterSet.cluster0);
    for (; !Kokkos::atomic_compare_exchange_strong(
                &_nucleationCrsClusterSets(id).cluster0,
                NetworkType::invalidIndex(), clusterSet.cluster0);
            ++id)
    {
    }
    _nucleationCrsClusterSets(id) = clusterSet;
}
}
}
}
}
