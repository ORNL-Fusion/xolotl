#pragma once

namespace xolotlCore
{
namespace experimental
{
namespace detail
{
template <typename TNetwork, typename TDerived>
class ReactionGenerator
{
public:
    struct Count
    {
    };

    struct Construct
    {
    };

    using Network = TNetwork;
    // using Types = detail::ReactionNetworkTypes<Network>;
    using ClusterData = typename Network::ClusterData;
    using ClusterDataRef = typename Network::ClusterDataRef;
    using Cluster = typename ClusterData::ClusterType;
    using ReactionType = typename Network::ReactionType;
    using ClusterSet = typename ReactionType::ClusterSet;
    using Subpaving = typename Network::Subpaving;

    ReactionGenerator(const Network& network)
        :
        _subpaving(network._subpaving),
        _clusterData(network._clusterData),
        _numDOFs(network.getDOF()),
        _clusterProdReactionCounts("Production Reaction Counts",
            _clusterData.numClusters),
        _clusterDissReactionCounts("Dissociation Reaction Counts",
            _clusterData.numClusters)
    {
    }

    void
    generateReactions()
    {
        auto numClusters = _clusterData.numClusters;
        auto diffusionFactor = _clusterData.diffusionFactor;
        auto generator = *static_cast<TDerived*>(this);
        using Range2D = Kokkos::MDRangePolicy<Kokkos::Rank<2>>;
        auto range2d = Range2D({0, 0}, {numClusters, numClusters});
        Kokkos::parallel_for(range2d,
                KOKKOS_LAMBDA (std::size_t i, std::size_t j) {
            if (j < i) {
                return;
            }
            if (diffusionFactor(i) == 0.0 && diffusionFactor(j) == 0.0) {
                return;
            }
            generator(i, j, Count{});
        });
        Kokkos::fence();

        setupCrs();
        setupReactionData();

        generator = *static_cast<TDerived*>(this);

        Kokkos::parallel_for(range2d,
                KOKKOS_LAMBDA (std::size_t i, std::size_t j) {
            if (j < i) {
                return;
            }
            if (diffusionFactor(i) == 0.0 && diffusionFactor(j) == 0.0) {
                return;
            }
            generator(i, j, Construct{});
        });
        Kokkos::fence();

        using RType = typename ReactionType::Type;
        auto reactionData = _reactionDataRef;
        auto clusterData = ClusterDataRef(_clusterData);
        auto prodReactions = _prodCrsReactions;
        auto prodClusterSets = _prodCrsClusterSets;
        Kokkos::parallel_for(_numProdReactions, KOKKOS_LAMBDA (std::size_t i) {
            prodReactions(i) = ReactionType(reactionData, clusterData, i,
                RType::production, prodClusterSets(i));
        });
        auto dissReactions = _dissCrsReactions;
        auto dissClusterSets = _dissCrsClusterSets;
        Kokkos::parallel_for(_numDissReactions, KOKKOS_LAMBDA (std::size_t i) {
            dissReactions(i) = ReactionType(reactionData, clusterData, i,
                RType::dissociation, dissClusterSets(i));
        });
        Kokkos::fence();
    }

    KOKKOS_INLINE_FUNCTION
    const Subpaving&
    getSubpaving() const
    {
        return _subpaving;
    }

    KOKKOS_INLINE_FUNCTION
    Cluster
    getCluster(std::size_t i) const
    {
        return _clusterData.getCluster(i);
    }

    void
    setupCrs()
    {
        _numProdReactions = Kokkos::get_crs_row_map_from_counts(
            _prodCrsRowMap, _clusterProdReactionCounts);

        _numDissReactions = Kokkos::get_crs_row_map_from_counts(
            _dissCrsRowMap, _clusterDissReactionCounts);

        _prodCrsClusterSets =
            Kokkos::View<ClusterSet*>("Production Cluster Sets",
                _numProdReactions);

        _dissCrsClusterSets =
            Kokkos::View<ClusterSet*>("Dissociation Cluster Sets",
                _numDissReactions);

        auto numReactions = _numProdReactions + _numDissReactions;
        _reactions = Kokkos::View<ReactionType*>(
                "Reactions", numReactions);
            // Kokkos::ViewAllocateWithoutInitializing("Reactions"), numReactions);

        _prodCrsReactions = Kokkos::subview(_reactions,
            std::make_pair((std::size_t)0, _numProdReactions));

        _dissCrsReactions = Kokkos::subview(_reactions,
            std::make_pair(_numProdReactions, numReactions));
    }

    void
    setupReactionData()
    {
        _reactionData = detail::ReactionData(_numProdReactions,
            _numDissReactions, _numDOFs, Network::getNumberOfSpeciesNoI(),
            _clusterData.gridSize);
        _reactionDataRef = detail::ReactionDataRef(_reactionData);
    }

    KOKKOS_INLINE_FUNCTION
    std::size_t
    getNumberOfClusters() const noexcept
    {
        return _clusterData.numClusters;
    }

    KOKKOS_INLINE_FUNCTION
    void
    addProductionReaction(Count, const ClusterSet& clusterSet) const
    {
        Kokkos::atomic_increment(
            &_clusterProdReactionCounts(clusterSet.cluster0));
    }

    KOKKOS_INLINE_FUNCTION
    void
    addProductionReaction(Construct, const ClusterSet& clusterSet) const
    {
        auto id = _prodCrsRowMap(clusterSet.cluster0);
        for (; !Kokkos::atomic_compare_exchange_strong(
                    &_prodCrsClusterSets(id).cluster0, Network::invalid,
                    clusterSet.cluster0);
                ++id)
        {
        }
        _prodCrsClusterSets(id) = clusterSet;
    }

    KOKKOS_INLINE_FUNCTION
    void
    addDissociationReaction(Count, const ClusterSet& clusterSet) const
    {
        Kokkos::atomic_increment(
            &_clusterDissReactionCounts(clusterSet.cluster1));
    }

    KOKKOS_INLINE_FUNCTION
    void
    addDissociationReaction(Construct, const ClusterSet& clusterSet) const
    {
        auto id = _dissCrsRowMap(clusterSet.cluster1);
        for (; !Kokkos::atomic_compare_exchange_strong(
                    &_dissCrsClusterSets(id).cluster0, Network::invalid,
                    clusterSet.cluster0);
                ++id)
        {
        }
        _dissCrsClusterSets(id) = clusterSet;
    }

    detail::ReactionData
    getReactionData() const
    {
        return _reactionData;
    }

    Kokkos::View<ReactionType*>
    getReactions() const
    {
        return _reactions;
    }

private:
    Subpaving _subpaving;
    ClusterData _clusterData;
    std::size_t _numDOFs;
    Kokkos::View<std::size_t*> _clusterProdReactionCounts;
    Kokkos::View<std::size_t*> _clusterDissReactionCounts;

    std::size_t _numProdReactions;
    std::size_t _numDissReactions;

    Kokkos::View<std::size_t*> _prodCrsRowMap;
    Kokkos::View<std::size_t*> _dissCrsRowMap;

    Kokkos::View<ClusterSet*> _prodCrsClusterSets;
    Kokkos::View<ClusterSet*> _dissCrsClusterSets;

    using ReactionSubView = decltype(
        Kokkos::subview(std::declval<Kokkos::View<ReactionType*>>(),
            Kokkos::ALL));
    ReactionSubView _prodCrsReactions;
    ReactionSubView _dissCrsReactions;

    detail::ReactionData _reactionData;
    detail::ReactionDataRef _reactionDataRef;
    Kokkos::View<ReactionType*> _reactions;
};
}
}
}
