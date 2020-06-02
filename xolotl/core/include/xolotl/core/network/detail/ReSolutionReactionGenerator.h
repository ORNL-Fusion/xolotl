#pragma once

#include <xolotl/core/network/detail/ReactionGenerator.h>

namespace xolotl
{
namespace core
{
namespace network
{
namespace detail
{
template <typename TBase>
class ReSolutionReactionGenerator : public TBase
{
public:
    using Superclass = TBase;
    using NetworkType = typename TBase::NetworkType;
    using NetworkTraits = ReactionNetworkTraits<NetworkType>;
    using ReSolutionReactionType =
        typename NetworkTraits::ReSolutionReactionType;
    using IndexType = typename NetworkType::IndexType;
    using IndexView = typename Superclass::IndexView;
    using ClusterSetSubView = typename Superclass::ClusterSetSubView;
    using Count = typename Superclass::Count;
    using Construct = typename Superclass::Construct;

    ReSolutionReactionGenerator(const NetworkType& network)
        :
        Superclass(network),
        _clusterReSoReactionCounts("ReSolution Reaction Counts",
            Superclass::getNumberOfClusters())
    {
    }

    IndexType
    getRowMapAndTotalReactionCount()
    {
        _numPrecedingReactions = Superclass::getRowMapAndTotalReactionCount();
        _numReSoReactions = Kokkos::get_crs_row_map_from_counts(_reSoCrsRowMap,
            _clusterReSoReactionCounts);

        _reSoReactions = Kokkos::View<ReSolutionReactionType*>(
            "ReSolution Reactions", _numReSoReactions);

        return _numPrecedingReactions + _numReSoReactions;
    }

    void
    setupCrsClusterSetSubView()
    {
        Superclass::setupCrsClusterSetSubView();
        _reSoCrsClusterSets = this->getClusterSetSubView(
            std::make_pair(_numPrecedingReactions,
                _numPrecedingReactions + _numReSoReactions));
    }

    KOKKOS_INLINE_FUNCTION
    void
    addReSolutionReaction(Count, const ClusterSet& clusterSet) const
    {
        if (!this->_clusterData.enableReSolution(0)) return;

        Kokkos::atomic_increment(
            &_clusterReSoReactionCounts(clusterSet.cluster0));
    }

    KOKKOS_INLINE_FUNCTION
    void
    addReSolutionReaction(Construct, const ClusterSet& clusterSet) const
    {
        if (!this->_clusterData.enableReSolution(0)) return;

        auto id = _reSoCrsRowMap(clusterSet.cluster0);
        for (; !Kokkos::atomic_compare_exchange_strong(
                    &_reSoCrsClusterSets(id).cluster0,
                    NetworkType::invalidIndex(), clusterSet.cluster0);
                ++id)
        {
        }
        _reSoCrsClusterSets(id) = clusterSet;
    }

    Kokkos::View<ReSolutionReactionType*>
    getReSolutionReactions() const
    {
        return _reSoReactions;
    }

    IndexType
    getNumberOfReSolutionReactions() const
    {
        return _reSoReactions.size();
    }

private:
    IndexView _clusterReSoReactionCounts;

    IndexType _numPrecedingReactions {};
    IndexType _numReSoReactions {};

    IndexView _reSoCrsRowMap;
    ClusterSetSubView _reSoCrsClusterSets;

    Kokkos::View<ReSolutionReactionType*> _reSoReactions;
};

template <typename TNetwork, typename TReaction, typename TBase>
struct WrapTypeSpecificReactionGenerator<TNetwork, TReaction, TBase,
    std::enable_if_t<
        std::is_base_of<ReSolutionReaction<TNetwork, TReaction>, TReaction>::value>>
{
    using Type = ReSolutionReactionGenerator<TBase>;
};
}
}
}
}
