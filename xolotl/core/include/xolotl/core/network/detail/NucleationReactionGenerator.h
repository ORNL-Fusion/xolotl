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
class NucleationReactionGenerator : public TBase
{
public:
    using Superclass = TBase;
    using NetworkType = typename TBase::NetworkType;
    using NetworkTraits = ReactionNetworkTraits<NetworkType>;
    using NucleationReactionType = typename NetworkTraits::NucleationReactionType;
    using IndexType = typename NetworkType::IndexType;
    using IndexView = typename Superclass::IndexView;
    using ClusterSetSubView = typename Superclass::ClusterSetSubView;
    using Count = typename Superclass::Count;
    using Construct = typename Superclass::Construct;

    NucleationReactionGenerator(const NetworkType& network)
        :
        Superclass(network),
        _clusterNucleationReactionCounts("Nucleation Reaction Counts",
            Superclass::getNumberOfClusters())
    {
    }

    IndexType
    getRowMapAndTotalReactionCount()
    {
        _numPrecedingReactions = Superclass::getRowMapAndTotalReactionCount();
        _numNucleationReactions = Kokkos::get_crs_row_map_from_counts(_nucleationCrsRowMap,
            _clusterNucleationReactionCounts);

        _nucleationReactions = Kokkos::View<NucleationReactionType*>("Nucleation Reactions",
            _numNucleationReactions);

        return _numPrecedingReactions + _numNucleationReactions;
    }

    void
    setupCrsClusterSetSubView()
    {
        Superclass::setupCrsClusterSetSubView();
        _nucleationCrsClusterSets = this->getClusterSetSubView(
            std::make_pair(_numPrecedingReactions,
                _numPrecedingReactions + _numNucleationReactions));
    }

    KOKKOS_INLINE_FUNCTION
    void
    addNucleationReaction(Count, const ClusterSet& clusterSet) const
    {
        Kokkos::atomic_increment(
            &_clusterNucleationReactionCounts(clusterSet.cluster0));
    }

    KOKKOS_INLINE_FUNCTION
    void
    addNucleationReaction(Construct, const ClusterSet& clusterSet) const
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

    Kokkos::View<NucleationReactionType*>
    getNucleationReactions() const
    {
        return _nucleationReactions;
    }

    IndexType
    getNumberOfNucleationReactions() const
    {
        return _nucleationReactions.size();
    }

private:
    IndexView _clusterNucleationReactionCounts;

    IndexType _numPrecedingReactions {};
    IndexType _numNucleationReactions {};

    IndexView _nucleationCrsRowMap;
    ClusterSetSubView _nucleationCrsClusterSets;

    Kokkos::View<NucleationReactionType*> _nucleationReactions;
};

template <typename TNetwork, typename TReaction, typename TBase>
struct WrapTypeSpecificReactionGenerator<TNetwork, TReaction, TBase,
    std::enable_if_t<
        std::is_base_of<NucleationReaction<TNetwork, TReaction>, TReaction>::value>>
{
    using Type = NucleationReactionGenerator<TBase>;
};
}
}
}
}
