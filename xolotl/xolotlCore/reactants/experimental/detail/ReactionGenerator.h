#pragma once

#include <utility>
#include <type_traits>

#include <experimental/ReactionNetworkTraits.h>

namespace xolotlCore
{
namespace experimental
{
namespace detail
{
template <typename TNetwork, typename TDerived>
class ReactionGeneratorBase
{
public:
    using NetworkType = TNetwork;
    using NetworkTraits = ReactionNetworkTraits<NetworkType>;
    using ClusterData = typename NetworkType::ClusterData;
    using ClusterDataRef = typename NetworkType::ClusterDataRef;
    using Cluster = typename ClusterData::ClusterType;
    using ProductionReactionType =
        typename NetworkTraits::ProductionReactionType;
    using DissociationReactionType =
        typename NetworkTraits::DissociationReactionType;
    using Subpaving = typename NetworkType::Subpaving;
    using IndexType = typename NetworkType::IndexType;
    using IndexView = Kokkos::View<IndexType*>;
    using ClusterSetView = Kokkos::View<ClusterSet*>;
    using ClusterSetSubView = decltype(Kokkos::subview(
        std::declval<ClusterSetView>(),
        std::declval<std::pair<IndexType, IndexType>>()));
    using Connectivity = typename NetworkType::Connectivity;

    struct Count
    {
    };

    struct Construct
    {
    };

    ReactionGeneratorBase(const TNetwork& network)
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
        auto generator = *(this->asDerived());
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

        generator = *(this->asDerived());

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

        auto reactionCollection = this->asDerived()->getReactionCollection();
        reactionCollection.construct(_reactionDataRef, _clusterData,
            _allClusterSets);

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

    IndexType
    getRowMapAndTotalReactionCount()
    {
        _numProdReactions = Kokkos::get_crs_row_map_from_counts(_prodCrsRowMap,
            _clusterProdReactionCounts);
        _numDissReactions = Kokkos::get_crs_row_map_from_counts(_dissCrsRowMap,
            _clusterDissReactionCounts);

        _prodReactions = Kokkos::View<ProductionReactionType*>(
            "Production Reactions", _numProdReactions);
        _dissReactions = Kokkos::View<DissociationReactionType*>(
            "Dissociation Reactions", _numDissReactions);

        return _numProdReactions + _numDissReactions;
    }

    ClusterSetSubView
    getClusterSetSubView(std::pair<IndexType, IndexType> indexRange)
    {
        return Kokkos::subview(_allClusterSets, indexRange);
    }

    void
    setupCrsClusterSetSubView()
    {
        _prodCrsClusterSets = getClusterSetSubView(
            std::make_pair(static_cast<IndexType>(0), _numProdReactions));

        _dissCrsClusterSets = getClusterSetSubView(
            std::make_pair(_numProdReactions,
                _numProdReactions + _numDissReactions));
    }

    void
    setupCrs()
    {
        auto numTotalReactions =
            this->asDerived()->getRowMapAndTotalReactionCount();
        _allClusterSets = ClusterSetView("Cluster Sets", numTotalReactions);
        this->asDerived()->setupCrsClusterSetSubView();
    }

    void
    setupReactionData()
    {
        _reactionData = detail::ReactionData(_numProdReactions,
            _numDissReactions, this->asDerived()->getNumberOfSinkReactions(),
            this->asDerived()->getNumberOfReSolutionReactions(),
            NetworkType::getNumberOfSpeciesNoI(), _clusterData.gridSize);
        _reactionDataRef = detail::ReactionDataRef(_reactionData);
    }

    IndexType
    getNumberOfSinkReactions() const noexcept
    {
        return 0;
    }

    IndexType
    getNumberOfReSolutionReactions() const noexcept
    {
        return 0;
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
        if (!_clusterData.enableStdReaction(0)) return;

        Kokkos::atomic_increment(
            &_clusterProdReactionCounts(clusterSet.cluster0));
    }

    KOKKOS_INLINE_FUNCTION
    void
    addProductionReaction(Construct, const ClusterSet& clusterSet) const
    {
        if (!_clusterData.enableStdReaction(0)) return;

        auto id = _prodCrsRowMap(clusterSet.cluster0);
        for (; !Kokkos::atomic_compare_exchange_strong(
                    &_prodCrsClusterSets(id).cluster0,
                    NetworkType::invalidIndex(), clusterSet.cluster0);
                ++id)
        {
        }
        _prodCrsClusterSets(id) = clusterSet;
    }

    KOKKOS_INLINE_FUNCTION
    void
    addDissociationReaction(Count, const ClusterSet& clusterSet) const
    {
        if (!_clusterData.enableStdReaction(0)) return;

        Kokkos::atomic_increment(
            &_clusterDissReactionCounts(clusterSet.cluster1));
    }

    KOKKOS_INLINE_FUNCTION
    void
    addDissociationReaction(Construct, const ClusterSet& clusterSet) const
    {
        if (!_clusterData.enableStdReaction(0)) return;

        auto id = _dissCrsRowMap(clusterSet.cluster1);
        for (; !Kokkos::atomic_compare_exchange_strong(
                    &_dissCrsClusterSets(id).cluster1, NetworkType::invalidIndex(),
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

    void
    generateConnectivity()
    {
        using RowMap = typename Connectivity::row_map_type;
        using Entries = typename Connectivity::entries_type;

        auto reactionCollection = this->asDerived()->getReactionCollection();

        Connectivity tmpConn;
        //Count connectivity entries
        //NOTE: We're using row_map for counts because
        //      Reaction::contributeConnectivity expects the connectivity CRS
        tmpConn.row_map =
            RowMap(Kokkos::ViewAllocateWithoutInitializing("tmp counts"),
                this->_numDOFs);
        // Even if there is no reaction each dof should connect with itself (for PETSc)
        Kokkos::parallel_for(this->_numDOFs, KOKKOS_LAMBDA (const IndexType i) {
            tmpConn.row_map(i) = 1;
        });
        reactionCollection.apply(DEVICE_LAMBDA (auto&& reaction) {
            reaction.contributeConnectivity(tmpConn);
        });

        Kokkos::fence();
        //Get row map
        auto counts = tmpConn.row_map;
        auto nEntries = Kokkos::get_crs_row_map_from_counts(tmpConn.row_map,
            counts);
        //Reset counts view
        counts = RowMap();
        //Initialize entries to invalid
        tmpConn.entries = Entries(
            Kokkos::ViewAllocateWithoutInitializing("connectivity entries"),
            nEntries);
        Kokkos::parallel_for(nEntries, KOKKOS_LAMBDA (IndexType i) {
            tmpConn.entries(i) = NetworkType::invalidIndex();
        });
        // Even if there is no reaction each dof should connect with itself (for PETSc)
        Kokkos::parallel_for(this->_numDOFs, KOKKOS_LAMBDA (const IndexType i) {
            auto id = tmpConn.row_map(i);
            for (; !Kokkos::atomic_compare_exchange_strong(
                        &tmpConn.entries(id), NetworkType::invalidIndex(), i);
                        ++id) {
                if (tmpConn.entries(id) == i) {
                    break;
                }
            }
        });
        //Fill entries (column ids)
        reactionCollection.apply(DEVICE_LAMBDA (auto&& reaction) {
            reaction.contributeConnectivity(tmpConn);
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
                    if (tmpConn.entries(j) == NetworkType::invalidIndex()) {
                        ret = j - jStart;
                        break;
                    }
                }
            }
            else {
                auto tmpStart = tmpConn.row_map(i);
                for (IndexType j = tmpStart; j < tmpConn.row_map(i+1); ++j) {
                    auto entry = tmpConn.entries(j);
                    if (entry == NetworkType::invalidIndex()) {
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

protected:
    TDerived*
    asDerived()
    {
        return static_cast<TDerived*>(this);
    }

protected:
    Subpaving _subpaving;
    ClusterData _clusterData;
    IndexType _numDOFs;
    IndexView _clusterProdReactionCounts;
    IndexView _clusterDissReactionCounts;

    IndexType _numProdReactions;
    IndexType _numDissReactions;

    Kokkos::View<IndexType*> _prodCrsRowMap;
    Kokkos::View<IndexType*> _dissCrsRowMap;

    ClusterSetView _allClusterSets;
    ClusterSetSubView _prodCrsClusterSets;
    ClusterSetSubView _dissCrsClusterSets;

    Kokkos::View<ProductionReactionType*> _prodReactions;
    Kokkos::View<DissociationReactionType*> _dissReactions;

    detail::ReactionData _reactionData;
    detail::ReactionDataRef _reactionDataRef;
};

template <typename TNetwork, typename TReaction,
    typename TReactionGeneratorParent, typename = void>
struct WrapTypeSpecificReactionGenerator
{
    // WrapTypeSpecificReactionGenerator()
    // {
    //     static_assert(false,
    //         "No type-specific reaction generator for this reaction type");
    // }
};

template <typename TReactionGeneratorParent, typename TExtraReactionTypes>
struct ReactionGeneratorTypeBuilderImpl;

template <typename TReactionGeneratorParent>
struct ReactionGeneratorTypeBuilderImpl<TReactionGeneratorParent, std::tuple<>>
{
    using NetworkType = typename TReactionGeneratorParent::NetworkType;
    using Type = TReactionGeneratorParent;
};

template <typename TReactionGeneratorParent, typename... TExtraReactions>
struct ReactionGeneratorTypeBuilderImpl<TReactionGeneratorParent,
    std::tuple<TExtraReactions...>>
{
    using NetworkType = typename TReactionGeneratorParent::NetworkType;
    using ExtraReactions = std::tuple<TExtraReactions...>;
    using FrontReaction = std::tuple_element_t<0, ExtraReactions>;
    using Type =
        typename WrapTypeSpecificReactionGenerator<NetworkType, FrontReaction,
            typename ReactionGeneratorTypeBuilderImpl<
                TReactionGeneratorParent, TuplePopFront<ExtraReactions>>::Type>
                    ::Type;
};

class NoneSuch { public: using NetworkType = void; };

template <typename TNetwork, typename TDerived>
struct ReactionGeneratorTypeBuilder
{
    using NetworkType = TNetwork;
    using ReactionTypes = ReactionTypeList<NetworkType>;
    using ReactionFirst = std::tuple_element_t<0, ReactionTypes>;
    using ReactionSecond = std::tuple_element_t<1, ReactionTypes>;

    static_assert(
        std::is_base_of<ProductionReaction<TNetwork, ReactionFirst>,
            ReactionFirst>::value,
        "First reaction type must be a ProductionReaction");

    static_assert(
        std::is_base_of<DissociationReaction<TNetwork, ReactionSecond>,
            ReactionSecond>::value,
        "Second reaction type must be a DissociationReaction");

    using ExtraReactionTypes = TuplePopFront<TuplePopFront<ReactionTypes>>;
    using Type = typename ReactionGeneratorTypeBuilderImpl<
        ReactionGeneratorBase<TNetwork, TDerived>,
        TupleReverse<ExtraReactionTypes>>::Type;
};

template <typename TNetwork, typename TDerived>
using ReactionGenerator =
    typename ReactionGeneratorTypeBuilder<TNetwork, TDerived>::Type;
}
}
}
