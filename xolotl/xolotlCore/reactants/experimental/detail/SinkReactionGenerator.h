#pragma once

#include <experimental/detail/ReactionGenerator.h>

namespace xolotlCore
{
namespace experimental
{
namespace detail
{
template <typename TBase>
class SinkReactionGenerator : public TBase
{
public:
    using Superclass = TBase;
    using NetworkType = typename TBase::NetworkType;
    using NetworkTraits = ReactionNetworkTraits<NetworkType>;
    using SinkReactionType = typename NetworkTraits::SinkReactionType;
    using IndexType = typename NetworkType::IndexType;
    using IndexView = typename Superclass::IndexView;
    using ClusterSetSubView = typename Superclass::ClusterSetSubView;
    using Count = typename Superclass::Count;
    using Construct = typename Superclass::Construct;

    SinkReactionGenerator(const NetworkType& network)
        :
        Superclass(network),
        _clusterSinkReactionCounts("Sink Reaction Counts",
            Superclass::getNumberOfClusters())
    {
    }

    IndexType
    getRowMapAndTotalReactionCount()
    {
        _numPrecedingReactions = Superclass::getRowMapAndTotalReactionCount();
        _numSinkReactions = Kokkos::get_crs_row_map_from_counts(_sinkCrsRowMap,
            _clusterSinkReactionCounts);

        _sinkReactions = Kokkos::View<SinkReactionType*>("Sink Reactions",
            _numSinkReactions);

        return _numPrecedingReactions + _numSinkReactions;
    }

    void
    setupCrsClusterSetSubView()
    {
        Superclass::setupCrsClusterSetSubView();
        _sinkCrsClusterSets = this->getClusterSetSubView(
            std::make_pair(_numPrecedingReactions,
                _numPrecedingReactions + _numSinkReactions));
    }

    KOKKOS_INLINE_FUNCTION
    void
    addSinkReaction(Count, const ClusterSet& clusterSet) const
    {
        if (!this->_clusterData.enableStdReaction(0)) return;

        Kokkos::atomic_increment(
            &_clusterSinkReactionCounts(clusterSet.cluster0));
    }

    KOKKOS_INLINE_FUNCTION
    void
    addSinkReaction(Construct, const ClusterSet& clusterSet) const
    {
        if (!this->_clusterData.enableStdReaction(0)) return;

        auto id = _sinkCrsRowMap(clusterSet.cluster0);
        for (; !Kokkos::atomic_compare_exchange_strong(
                    &_sinkCrsClusterSets(id).cluster0,
                    NetworkType::invalidIndex(), clusterSet.cluster0);
                ++id)
        {
        }
        _sinkCrsClusterSets(id) = clusterSet;
    }

    Kokkos::View<SinkReactionType*>
    getSinkReactions() const
    {
        return _sinkReactions;
    }

    IndexType
    getNumberOfSinkReactions() const
    {
        return _sinkReactions.size();
    }

private:
    IndexView _clusterSinkReactionCounts;

    IndexType _numPrecedingReactions {};
    IndexType _numSinkReactions {};

    IndexView _sinkCrsRowMap;
    ClusterSetSubView _sinkCrsClusterSets;

    Kokkos::View<SinkReactionType*> _sinkReactions;
};

template <typename TNetwork, typename TReaction, typename TBase>
struct WrapTypeSpecificReactionGenerator<TNetwork, TReaction, TBase,
    std::enable_if_t<
        std::is_base_of<SinkReaction<TNetwork, TReaction>, TReaction>::value>>
{
    using Type = SinkReactionGenerator<TBase>;
};
}
}
}
