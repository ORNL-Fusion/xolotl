#pragma once

namespace xolotlCore
{
namespace experimental
{
namespace detail
{
template <typename TNetwork, typename TDerived>
class ReactionGeneratorRange1D
{
public:
    struct Count
    {
    };

    struct Construct
    {
    };

    using Network = TNetwork;
    using ClusterData = typename Network::ClusterData;
    using Cluster = typename ClusterData::ClusterType;
    using ProductionReactionType = typename Network::ProductionReactionType;
    using DissociationReactionType = typename Network::DissociationReactionType;
    using Subpaving = typename Network::Subpaving;
    using IndexType = typename Network::IndexType;

    ReactionGeneratorRange1D(const Subpaving& subpaving,
            const ClusterData& clusterData, IndexType numDOFs)
        :
        _subpaving(subpaving),
        _clusterData(clusterData),
        _numDOFs(numDOFs),
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
        Kokkos::parallel_for(numClusters, KOKKOS_LAMBDA (IndexType j) {
            for (IndexType i = 0; i <= j; ++i) {
                if (diffusionFactor(i) == 0.0 && diffusionFactor(j) == 0.0) {
                    continue;
                }
                generator(i, j, Count{});
            }
        });
        Kokkos::fence();

        setupCrs();
        setupReactionData();

        generator = *static_cast<TDerived*>(this);

        Kokkos::parallel_for(numClusters, KOKKOS_LAMBDA (IndexType j) {
            for (IndexType i = 0; i <= j; ++i) {
                if (diffusionFactor(i) == 0.0 && diffusionFactor(j) == 0.0) {
                    continue;
                }
                generator(i, j, Construct{});
            }
        });
        Kokkos::fence();

        static_cast<TDerived*>(this)->generateConnectivity();
    }

    KOKKOS_INLINE_FUNCTION
    const Subpaving&
    getSubpaving() const
    {
        return _subpaving;
    }

    KOKKOS_INLINE_FUNCTION
    Cluster
    getCluster(IndexType i) const
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

        _prodReactions = Kokkos::View<ProductionReactionType*>(
            "Production Reactions", _numProdReactions);
        _dissReactions = Kokkos::View<DissociationReactionType*>(
            "Dissociation Reactions", _numDissReactions);
    }

    void
    setupReactionData()
    {
        _reactionData = detail::ReactionData(_numProdReactions,
            _numDissReactions, Network::getNumberOfSpeciesNoI(),
            _clusterData.gridSize);
        _reactionDataRef = detail::ReactionDataRef(_reactionData);
    }

    KOKKOS_INLINE_FUNCTION
    IndexType
    getNumberOfClusters() const noexcept
    {
        return _clusterData.numClusters;
    }

    KOKKOS_INLINE_FUNCTION
    void
    addProductionReaction(Count, const ClusterSet& clusterSet) const
    {
        ++_clusterProdReactionCounts(clusterSet.cluster1);
    }

    KOKKOS_INLINE_FUNCTION
    void
    addProductionReaction(Construct, const ClusterSet& clusterSet) const
    {
        auto id = _prodCrsRowMap(clusterSet.cluster1);
        _prodReactions(id) = ProductionReactionType(_reactionDataRef,
            _clusterData, id, clusterSet);
        ++_prodCrsRowMap(clusterSet.cluster1);
    }

    KOKKOS_INLINE_FUNCTION
    void
    addDissociationReaction(Count, const ClusterSet& clusterSet) const
    {
        ++_clusterDissReactionCounts(clusterSet.cluster2);
    }

    KOKKOS_INLINE_FUNCTION
    void
    addDissociationReaction(Construct, const ClusterSet& clusterSet) const
    {
        auto id = _dissCrsRowMap(clusterSet.cluster2);
        _dissReactions(id) = DissociationReactionType(_reactionDataRef,
            _clusterData, id + _numProdReactions, clusterSet);
        ++_dissCrsRowMap(clusterSet.cluster2);
    }

    detail::ReactionData
    getReactionData() const
    {
        return _reactionData;
    }

    Kokkos::View<ProductionReactionType*>
    getProductionReactions() const
    {
        return _prodReactions;
    }

    Kokkos::View<DissociationReactionType*>
    getDissociationReactions() const
    {
        return _dissReactions;
    }

protected:
    Subpaving _subpaving;
    ClusterData _clusterData;
    IndexType _numDOFs;
    Kokkos::View<IndexType*> _clusterProdReactionCounts;
    Kokkos::View<IndexType*> _clusterDissReactionCounts;

    IndexType _numProdReactions;
    IndexType _numDissReactions;

    Kokkos::View<IndexType*> _prodCrsRowMap;
    Kokkos::View<IndexType*> _dissCrsRowMap;

    Kokkos::View<ProductionReactionType*> _prodReactions;
    Kokkos::View<DissociationReactionType*> _dissReactions;

    detail::ReactionData _reactionData;
    detail::ReactionDataRef _reactionDataRef;
};

template <typename TNetwork, typename TDerived>
class ReactionGeneratorRange2D
{
public:
    struct Count
    {
    };

    struct Construct
    {
    };

    using Network = TNetwork;
    using ClusterData = typename Network::ClusterData;
    using ClusterDataRef = typename Network::ClusterDataRef;
    using Cluster = typename ClusterData::ClusterType;
    using ProductionReactionType = typename Network::ProductionReactionType;
    using DissociationReactionType = typename Network::DissociationReactionType;
    using Subpaving = typename Network::Subpaving;
    using IndexType = typename Network::IndexType;

    ReactionGeneratorRange2D(const Subpaving& subpaving,
            const ClusterData& clusterData, IndexType numDOFs)
        :
        _subpaving(subpaving),
        _clusterData(clusterData),
        _numDOFs(numDOFs),
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
                KOKKOS_LAMBDA (IndexType i, IndexType j) {
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
                KOKKOS_LAMBDA (IndexType i, IndexType j) {
            if (j < i) {
                return;
            }
            if (diffusionFactor(i) == 0.0 && diffusionFactor(j) == 0.0) {
                return;
            }
            generator(i, j, Construct{});
        });
        Kokkos::fence();

        auto reactionData = _reactionDataRef;
        auto clusterData = ClusterDataRef(_clusterData);
        auto prodReactions = _prodReactions;
        auto prodClusterSets = _prodCrsClusterSets;
        Kokkos::parallel_for(_numProdReactions, KOKKOS_LAMBDA (IndexType i) {
            prodReactions(i) = ProductionReactionType(reactionData, clusterData,
                i, prodClusterSets(i));
        });
        auto dissReactions = _dissReactions;
        auto dissClusterSets = _dissCrsClusterSets;
        auto nProdReactions = _numProdReactions;
        Kokkos::parallel_for(_numDissReactions, KOKKOS_LAMBDA (IndexType i) {
            dissReactions(i) = DissociationReactionType(reactionData,
                clusterData, i + nProdReactions, dissClusterSets(i));
        });
        Kokkos::fence();

        static_cast<TDerived*>(this)->generateConnectivity();
    }

    KOKKOS_INLINE_FUNCTION
    const Subpaving&
    getSubpaving() const
    {
        return _subpaving;
    }

    KOKKOS_INLINE_FUNCTION
    Cluster
    getCluster(IndexType i) const
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

        _prodReactions = Kokkos::View<ProductionReactionType*>(
            "Production Reactions", _numProdReactions);
        _dissReactions = Kokkos::View<DissociationReactionType*>(
            "Dissociation Reactions", _numDissReactions);
    }

    void
    setupReactionData()
    {
        _reactionData = detail::ReactionData(_numProdReactions,
            _numDissReactions, Network::getNumberOfSpeciesNoI(),
            _clusterData.gridSize);
        _reactionDataRef = detail::ReactionDataRef(_reactionData);
    }

    KOKKOS_INLINE_FUNCTION
    IndexType
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
                    &_prodCrsClusterSets(id).cluster0, Network::invalidIndex(),
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
                    &_dissCrsClusterSets(id).cluster1, Network::invalidIndex(),
                    clusterSet.cluster1);
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

    Kokkos::View<ProductionReactionType*>
    getProductionReactions() const
    {
        return _prodReactions;
    }

    Kokkos::View<DissociationReactionType*>
    getDissociationReactions() const
    {
        return _dissReactions;
    }

protected:
    Subpaving _subpaving;
    ClusterData _clusterData;
    IndexType _numDOFs;
    Kokkos::View<IndexType*> _clusterProdReactionCounts;
    Kokkos::View<IndexType*> _clusterDissReactionCounts;

    IndexType _numProdReactions;
    IndexType _numDissReactions;

    Kokkos::View<IndexType*> _prodCrsRowMap;
    Kokkos::View<IndexType*> _dissCrsRowMap;

    Kokkos::View<ClusterSet*> _prodCrsClusterSets;
    Kokkos::View<ClusterSet*> _dissCrsClusterSets;

    Kokkos::View<ProductionReactionType*> _prodReactions;
    Kokkos::View<DissociationReactionType*> _dissReactions;

    detail::ReactionData _reactionData;
    detail::ReactionDataRef _reactionDataRef;
};

template <typename TNetwork, typename TDerived>
using ReactionGeneratorImpl =
#if 1
ReactionGeneratorRange2D<TNetwork, TDerived>;
#else
ReactionGeneratorRange1D<TNetwork, TDerived>;
#endif

template <typename TNetwork, typename TDerived>
class ReactionGenerator : public ReactionGeneratorImpl<TNetwork, TDerived>
{
public:
    using Superclass = ReactionGeneratorImpl<TNetwork, TDerived>;
    using Network = typename Superclass::Network;
    using Connectivity = typename Network::Connectivity;
    using IndexType = typename Superclass::IndexType;

    using Superclass::Superclass;

    ReactionGenerator(const TNetwork& network)
        :
        Superclass(network._subpaving, network._clusterData, network.getDOF())
    {
    }

    void
    generateConnectivity()
    {
        using RowMap = typename Connectivity::row_map_type;
        using Entries = typename Connectivity::entries_type;

        auto prodReactions = this->_prodReactions;
        auto dissReactions = this->_dissReactions;
        auto numProdReactions = prodReactions.size();
        auto numReactions = numProdReactions + dissReactions.size();
        Connectivity tmpConn;
        //Count connectivity entries
        //NOTE: We're using row_map for counts because
        //      Reaction::contributeConnectivity expects the connectivity CRS
        tmpConn.row_map = RowMap("tmp counts", this->_numDOFs);
        // Even if there is no reaction each dof should connect with itself (for PETSc)
        Kokkos::parallel_for(this->_numDOFs, KOKKOS_LAMBDA (const IndexType i) {
            Kokkos::atomic_increment(&tmpConn.row_map(i));
        });
        Kokkos::parallel_for(numReactions, KOKKOS_LAMBDA (IndexType i) {
            if (i < numProdReactions) {
                prodReactions(i).contributeConnectivity(tmpConn);
            }
            else {
                dissReactions(i-numProdReactions).contributeConnectivity(tmpConn);
            }
        });
        Kokkos::fence();
        //Get row map
        auto counts = tmpConn.row_map;
        auto nEntries =
            Kokkos::get_crs_row_map_from_counts(tmpConn.row_map, counts);
        //Reset counts view
        counts = RowMap();
        //Initialize entries to invalid
        tmpConn.entries = Entries(
            Kokkos::ViewAllocateWithoutInitializing("connectivity entries"),
            nEntries);
        Kokkos::parallel_for(nEntries, KOKKOS_LAMBDA (IndexType i) {
            tmpConn.entries(i) = Network::invalidIndex();
        });
        // Even if there is no reaction each dof should connect with itself (for PETSc)
        Kokkos::parallel_for(this->_numDOFs, KOKKOS_LAMBDA (const IndexType i) {
            auto id = tmpConn.row_map(i);
            for (; !Kokkos::atomic_compare_exchange_strong(
                        &tmpConn.entries(id), Network::invalidIndex(), i); ++id) {
                if (tmpConn.entries(id) == i) {
                    break;
                }
            }
        });
        //Fill entries (column ids)
        Kokkos::parallel_for(numReactions, KOKKOS_LAMBDA (IndexType i) {
            if (i < numProdReactions) {
                prodReactions(i).contributeConnectivity(tmpConn);
            }
            else {
                dissReactions(i-numProdReactions).contributeConnectivity(tmpConn);
            }
        });
        Kokkos::fence();

        //Shrink to fit
        Connectivity connectivity;
        Kokkos::count_and_fill_crs(connectivity, this->_numDOFs,
                KOKKOS_LAMBDA (IndexType i, IndexType* fill) {
            IndexType ret = 0;
            if (fill == nullptr) {
                auto jStart = tmpConn.row_map(i);
                auto jEnd = tmpConn.row_map(i+1);
                ret = jEnd - jStart;
                for (IndexType j = jStart; j < jEnd; ++j) {
                    if (tmpConn.entries(j) == Network::invalidIndex()) {
                        ret = j - jStart;
                        break;
                    }
                }
            }
            else {
                auto tmpStart = tmpConn.row_map(i);
                for (IndexType j = tmpStart; j < tmpConn.row_map(i+1); ++j) {
                    auto entry = tmpConn.entries(j);
                    if (entry == Network::invalidIndex()) {
                        break;
                    }
                    fill[j - tmpStart] = entry;
                }
            }
            return ret;
        });
        nEntries = connectivity.entries.extent(0);

        this->_reactionData.connectivity = connectivity;
    }
};
}
}
}
