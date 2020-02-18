#pragma once

#include <cstddef>
#include <cstdint>
#include <type_traits>

#include <Kokkos_Core.hpp>
#include <Kokkos_Crs.hpp>

#include <plsm/Subpaving.h>
#include <plsm/refine/RegionDetector.h>

#include <IOptions.h>
#include <Options.h>

#include <experimental/Cluster.h>
#include <experimental/IReactionNetwork.h>
#include <experimental/Reaction.h>
#include <experimental/SpeciesEnumSequence.h>

namespace xolotlCore
{
namespace experimental
{
namespace detail
{
template <typename TImpl>
struct ReactionNetworkWorker;

template <typename TImpl, typename TDerived>
class ReactionGenerator;

template <typename TData>
class UpperTriangle
{
public:
    UpperTriangle(const std::string& label, std::size_t N)
        :
        _N(N),
        _data(label, ((N+1)*N)/2)
    {
    }

    KOKKOS_INLINE_FUNCTION
    std::size_t
    size() const noexcept
    {
        return _data.extent(0);
    }

    KOKKOS_INLINE_FUNCTION
    TData&
    operator()(std::size_t i, std::size_t j) const
    {
        assert(j >= i);
        assert(i < _N);
        assert(j < _N);
        auto id = i*_N + j - ((i+1)*i)/2;
        return _data(id);
    }

    KOKKOS_INLINE_FUNCTION
    TData&
    operator()(std::size_t i) const
    {
        assert(i < size());
        return _data(i);
    }

private:
    std::size_t _N;
    Kokkos::View<TData*> _data;
};
}

template <typename TImpl>
class ReactionNetwork : public IReactionNetwork
{
    friend class detail::ReactionNetworkWorker<TImpl>;
    template <typename, typename>
    friend class detail::ReactionGenerator;

    static constexpr auto invalid = plsm::invalid<std::size_t>;

public:
    using Traits = ReactionNetworkTraits<TImpl>;
    using Species = typename Traits::Species;

private:
    using Types = detail::ReactionNetworkTypes<TImpl>;

    static constexpr std::size_t numSpecies = Traits::numSpecies;

public:
    using SpeciesSequence = SpeciesEnumSequence<Species, numSpecies>;
    using SpeciesRange = EnumSequenceRange<Species, numSpecies>;
    using ReactionType = typename Traits::ReactionType;
    using ClusterGenerator = typename Traits::ClusterGenerator;
    using AmountType = typename IReactionNetwork::AmountType;
    using Subpaving = typename Types::Subpaving;
    using SubdivisionRatio = plsm::SubdivisionRatio<numSpecies>;
    using Composition = typename Types::Composition;
    using Region = typename Types::Region;
    using Ival = typename Region::IntervalType;
    using ConcentrationsView = typename IReactionNetwork::ConcentrationsView;
    using FluxesView = typename IReactionNetwork::FluxesView;
    using ConnectivityView = typename IReactionNetwork::ConnectivityView;
    using SparseFillMap = typename IReactionNetwork::SparseFillMap;
    using ClusterData = typename Types::ClusterData;
    using ClusterDataMirror = typename Types::ClusterDataMirror;
    using ClusterDataRef = typename Types::ClusterDataRef;
    using InverseMap = Kokkos::View<std::size_t**>;

    template <typename PlsmContext>
    using Cluster = Cluster<TImpl, PlsmContext>;

    ReactionNetwork() = default;

    ReactionNetwork(const Subpaving& subpaving, std::size_t gridSize,
        const IOptions& options);

    ReactionNetwork(const Subpaving& subpaving, std::size_t gridSize);

    ReactionNetwork(const std::vector<AmountType>& maxSpeciesAmounts,
        const std::vector<SubdivisionRatio>& subdivisionRatios,
        std::size_t gridSize, const IOptions& options);

    ReactionNetwork(const std::vector<AmountType>& maxSpeciesAmounts,
        std::size_t gridSize, const IOptions& options);

    KOKKOS_INLINE_FUNCTION
    static
    constexpr std::size_t
    getNumberOfSpecies() noexcept
    {
        return SpeciesSequence::size();
    }

    KOKKOS_INLINE_FUNCTION
    static
    constexpr std::size_t
    getNumberOfSpeciesNoI() noexcept
    {
        return SpeciesSequence::sizeNoI();
    }

    KOKKOS_INLINE_FUNCTION
    static
    constexpr SpeciesRange
    getSpeciesRange() noexcept
    {
        return SpeciesRange{};
    }

    KOKKOS_INLINE_FUNCTION
    static
    constexpr SpeciesRange
    getSpeciesRangeNoI() noexcept
    {
        return SpeciesRange(SpeciesSequence::first(),
            SpeciesSequence::lastNoI());
    }

    void
    setLatticeParameter(double latticeParameter) override;

    void
    setImpurityRadius(double impurityRadius) noexcept override
    {
        this->_impurityRadius =
            asDerived()->checkImpurityRadius(impurityRadius);
    }

    void
    setTemperatures(const std::vector<double>& gridTemperatures) override;

    void
    syncClusterDataOnHost() override
    {
        _subpaving.syncTiles(plsm::onHost);
        auto mirror = ClusterDataMirror(_subpaving, this->_gridSize);
        Kokkos::deep_copy(mirror.atomicVolume, _clusterData.atomicVolume);
        Kokkos::deep_copy(mirror.temperature, _clusterData.temperature);
        Kokkos::deep_copy(mirror.momentIds, _clusterData.momentIds);
        Kokkos::deep_copy(mirror.reactionRadius, _clusterData.reactionRadius);
        Kokkos::deep_copy(mirror.formationEnergy, _clusterData.formationEnergy);
        Kokkos::deep_copy(mirror.migrationEnergy, _clusterData.migrationEnergy);
        Kokkos::deep_copy(mirror.diffusionFactor, _clusterData.diffusionFactor);
        Kokkos::deep_copy(mirror.diffusionCoefficient, _clusterData.diffusionCoefficient);
        _clusterDataMirror = mirror;
    }

    KOKKOS_INLINE_FUNCTION
    Cluster<plsm::OnDevice>
    findCluster(const Composition& comp, plsm::OnDevice context)
    {
        return Cluster<plsm::OnDevice>(_clusterData,
            _subpaving.findTileId(comp, context));
    }

    Cluster<plsm::OnHost>
    findCluster(const Composition& comp, plsm::OnHost context)
    {
        return Cluster<plsm::OnHost>(_clusterDataMirror,
            _subpaving.findTileId(comp, context));
    }

    KOKKOS_INLINE_FUNCTION
    auto
    findCluster(const Composition& comp)
    {
        return findCluster(comp, plsm::onDevice);
    }

    ClusterCommon<plsm::OnHost>
    getClusterCommon(std::size_t clusterId) const override
    {
        return ClusterCommon<plsm::OnHost>(_clusterDataMirror, clusterId);
    }

    KOKKOS_INLINE_FUNCTION
    Cluster<plsm::OnDevice>
    getCluster(std::size_t clusterId, plsm::OnDevice)
    {
        return Cluster<plsm::OnDevice>(_clusterData, clusterId);
    }

    Cluster<plsm::OnHost>
    getCluster(std::size_t clusterId, plsm::OnHost)
    {
        return Cluster<plsm::OnHost>(_clusterDataMirror, clusterId);
    }

    KOKKOS_INLINE_FUNCTION
    auto
    getCluster(std::size_t clusterId)
    {
        return getCluster(clusterId, plsm::onDevice);
    }

    KOKKOS_INLINE_FUNCTION
    Subpaving&
    getSubpaving()
    {
        return _subpaving;
    }

    void
    computeAllFluxes(ConcentrationsView concentrations, FluxesView fluxes,
        std::size_t gridIndex) override;

    void
    computeAllPartials(ConcentrationsView concentrations,
        Kokkos::View<double*> values, std::size_t gridIndex) override;

    double
    getLargestRate() override;

    std::size_t
    getDiagonalFill(SparseFillMap& fillMap) override;

    /**
     * Get the total concentration of a given type of clusters.
     *
     * @param concentration The vector of concentrations
     * @param type The type of atom we want the concentration of
     * @param minSize The minimum number of atom to start counting
     * @return The total accumulated concentration times the size of the cluster
     */
    double
    getTotalAtomConcentration(ConcentrationsView concentrations,
            Species type, std::size_t minSize = 0);

    /**
     * Get the total concentration of a given type of clusters only if it is trapped in a vacancy.
     *
     * @param concentration The vector of concentrations
     * @param type The type of atom we want the concentration of
     * @param minSize The minimum number of atom to start counting
     * @return The total accumulated concentration times the size of the cluster
     */
    double
    getTotalTrappedAtomConcentration(ConcentrationsView concentrations,
            Species type, std::size_t minSize = 0);

private:
    KOKKOS_INLINE_FUNCTION
    TImpl*
    asDerived()
    {
        return static_cast<TImpl*>(this);
    }

    void
    defineMomentIds();

    void
    generateClusterData(const ClusterGenerator& generator);

    void
    defineReactions();

    void
    updateDiffusionCoefficients();

    KOKKOS_INLINE_FUNCTION
    double
    getTemperature(std::size_t gridIndex) const noexcept
    {
        return _clusterData.temperature(gridIndex);
    }

    KOKKOS_INLINE_FUNCTION
    std::size_t
    getReactionIndex(std::size_t rowId, std::size_t colId) const noexcept
    {
        return _reactionData.inverseMap(rowId, colId);
    }

private:
    Subpaving _subpaving;

    ClusterData _clusterData;
    ClusterDataMirror _clusterDataMirror;

    detail::ReactionData _reactionData;
    Kokkos::View<ReactionType*> _reactions;

    detail::ReactionNetworkWorker<TImpl> _worker;
};


namespace detail
{
template <typename TImpl>
struct ReactionNetworkWorker
{
    class ExclusiveScanFunctor
    {
    public:
        ExclusiveScanFunctor(Kokkos::View<std::size_t*> data)
            :
            _data(data)
        {
        }

        KOKKOS_INLINE_FUNCTION
        void
        operator()(std::size_t index, std::size_t& update, const bool finalPass)
            const
        {
            const auto temp = _data(index);
            if (finalPass) {
                _data(index) = update;
            }
            update += temp;
        }

        KOKKOS_INLINE_FUNCTION
        void
        init(std::size_t& update) const
        {
          update = 0;
        }

        KOKKOS_INLINE_FUNCTION
        void
        join(volatile std::size_t& update, std::size_t input) const
        {
          update += input;
        }

    private:
        Kokkos::View<std::size_t*> _data;
    };

    using Network = ReactionNetwork<TImpl>;
    using Types = ReactionNetworkTypes<TImpl>;
    using ClusterData = typename Types::ClusterData;
    using ClusterDataRef = typename Types::ClusterDataRef;

    Network& _nw;

    ReactionNetworkWorker(Network& network)
        :
        _nw(network)
    {
    }

    void
    updateDiffusionCoefficients();

    void
    generateClusterData(const typename Network::ClusterGenerator& generator);

    void
    defineMomentIds();

    struct ClusterSetsViewPair
    {
        using ClusterSet = typename Network::ReactionType::ClusterSet;
        Kokkos::View<ClusterSet*> prodClusterSets;
        Kokkos::View<ClusterSet*> dissClusterSets;
    };

    void
    defineReactions();

    std::size_t
    getDiagonalFill(typename Network::SparseFillMap& fillMap);
};

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
    using Types = detail::ReactionNetworkTypes<Network>;
    using ClusterData = typename Types::ClusterData;
    using Cluster = typename ClusterData::ClusterType;
    using ReactionType = typename Network::ReactionType;
    using ClusterSet = typename ReactionType::ClusterSet;

    ReactionGenerator(const Network& network)
        :
        _clusterData(network._clusterData),
        _numDOFs(network.getDOF()),
        _clusterProdReactionCounts("Production Reaction Counts",
            _clusterData.numClusters),
        _clusterDissReactionCounts("Dissociation Reaction Counts",
            _clusterData.numClusters)
    {
    }

    KOKKOS_INLINE_FUNCTION
    Count
    countTag() const noexcept
    {
        return {};
    }

    KOKKOS_INLINE_FUNCTION
    Construct
    constructTag() const noexcept
    {
        return {};
    }

    void
    generateReactions()
    {
        auto numClusters = _clusterData.numClusters;
        auto diffusionFactor = _clusterData.diffusionFactor;
        auto generator = *static_cast<TDerived*>(this);
        Kokkos::parallel_for(numClusters, KOKKOS_LAMBDA (std::size_t i) {
            std::size_t prodCount {};
            std::size_t dissCount {};
            for (std::size_t j = i; j < numClusters; ++j) {
                if (diffusionFactor(i) == 0.0 && diffusionFactor(j) == 0.0) {
                    continue;
                }
                generator(i, j, prodCount, dissCount, generator.countTag());
            }
        });
        Kokkos::fence();

        setupCrs();
        setupReactionData();

        generator = *static_cast<TDerived*>(this);

        Kokkos::parallel_for(numClusters, KOKKOS_LAMBDA (std::size_t i) {
            std::size_t prodCount {};
            std::size_t dissCount {};
            for (std::size_t j = i; j < numClusters; ++j) {
                if (diffusionFactor(i) == 0.0 && diffusionFactor(j) == 0.0) {
                    continue;
                }
                generator(i, j, prodCount, dissCount, generator.constructTag());
            }
        });
        Kokkos::fence();
    }

    KOKKOS_INLINE_FUNCTION
    std::size_t
    getNumberOfClusters() const noexcept
    {
        return _clusterData.numClusters;
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

        auto numReactions = _numProdReactions + _numDissReactions;
        _reactions = Kokkos::View<ReactionType*>("Reactions", numReactions);

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
    void
    addProductionReaction(Count, const ClusterSet& clusterSet,
        std::size_t& count) const
    {
        ++_clusterProdReactionCounts(clusterSet.cluster0);
    }

    KOKKOS_INLINE_FUNCTION
    void
    addProductionReaction(Construct, const ClusterSet& clusterSet,
        std::size_t& count) const
    {
        using RType = typename ReactionType::Type;
        auto id = _prodCrsRowMap(clusterSet.cluster0) + count;
        _prodCrsReactions(id) = ReactionType(_reactionDataRef, _clusterData, id,
            RType::production, clusterSet);
        ++count;
    }

    KOKKOS_INLINE_FUNCTION
    void
    addDissociationReaction(Count, const ClusterSet& clusterSet,
        std::size_t& count) const
    {
        ++_clusterDissReactionCounts(clusterSet.cluster1);
    }

    KOKKOS_INLINE_FUNCTION
    void
    addDissociationReaction(Construct, const ClusterSet& clusterSet,
        std::size_t& count) const
    {
        using RType = typename ReactionType::Type;
        auto id = _dissCrsRowMap(clusterSet.cluster1) + count;
        _dissCrsReactions(id) = ReactionType(_reactionDataRef, _clusterData,
            id + _numProdReactions, RType::dissociation, clusterSet);
        ++count;
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
    ClusterData _clusterData;
    std::size_t _numDOFs;
    Kokkos::View<std::size_t*> _clusterProdReactionCounts;
    Kokkos::View<std::size_t*> _clusterDissReactionCounts;

    std::size_t _numProdReactions;
    std::size_t _numDissReactions;

    Kokkos::View<std::size_t*> _prodCrsRowMap;
    Kokkos::View<std::size_t*> _dissCrsRowMap;
    using ReactionSubView = decltype(
        Kokkos::subview(std::declval<Kokkos::View<ReactionType*>>(),
            Kokkos::ALL));
    ReactionSubView _prodCrsReactions;
    ReactionSubView _dissCrsReactions;
    // Kokkos::Crs<ClusterSet, detail::DefaultMemorySpace> _prodCrs;
    // Kokkos::Crs<ClusterSet, detail::DefaultMemorySpace> _dissCrs;

    detail::ReactionData _reactionData;
    detail::ReactionDataRef _reactionDataRef;
    Kokkos::View<ReactionType*> _reactions;
};
}
}
}

#include <experimental/ReactionNetwork.inl>
