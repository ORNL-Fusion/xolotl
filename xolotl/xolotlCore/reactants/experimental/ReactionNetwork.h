#pragma once

#include <cstddef>
#include <cstdint>
#include <type_traits>

#include <plsm/Subpaving.h>
#include <plsm/refine/RegionDetector.h>

#include <IOptions.h>
#include <Options.h>

#include <experimental/EnumSequence.h>
#include <experimental/Cluster.h>
#include <experimental/Reaction.h>

namespace xolotlCore
{
namespace experimental
{
namespace detail
{
template <typename TImpl>
struct ReactionNetworkWorker;
}

template <typename TImpl>
class ReactionNetwork
{
    friend class detail::ReactionNetworkWorker<TImpl>;

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
    using AmountType = typename Types::AmountType;
    using Subpaving = typename Types::Subpaving;
    using SubdivisionRatio = typename Subpaving::SubdivisionRatioType;
    using Composition = typename Types::Composition;
    using Region = typename Types::Region;
    using Ival = typename Region::IntervalType;
    using ConcentrationsView = Kokkos::View<double*, Kokkos::MemoryUnmanaged>;
    using FluxesView = Kokkos::View<double*, Kokkos::MemoryUnmanaged>;
    using ConnectivityView = typename Types::ConnectivityView;
    using SparseFillMap = std::unordered_map<int, std::vector<int>>;
    using ClusterData = typename Types::ClusterData;
    using ClusterDataMirror = typename Types::ClusterDataMirror;
    using ClusterDataRef = typename Types::ClusterDataRef;
    using InverseMap = Kokkos::View<std::size_t**>;

    template <typename PlsmContext>
    using Cluster = Cluster<Subpaving, PlsmContext>;

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

    KOKKOS_INLINE_FUNCTION
    std::size_t
    getDOF() const noexcept
    {
        return _numDOFs;
    }

    KOKKOS_INLINE_FUNCTION
    double
    getLatticeParameter() const noexcept
    {
        return _latticeParameter;
    }

    KOKKOS_INLINE_FUNCTION
    double
    getAtomicVolume() const noexcept
    {
        return _atomicVolume;
    }

    void
    setLatticeParameter(double latticeParameter);

    KOKKOS_INLINE_FUNCTION
    double
    getInterstitialBias() const noexcept
    {
        return _interstitialBias;
    }

    void
    setInterstitialBias(double interstitialBias) noexcept
    {
        _interstitialBias = interstitialBias;
    }

    KOKKOS_INLINE_FUNCTION
    double
    getImpurityRadius() const noexcept
    {
        return _impurityRadius;
    }

    void
    setImpurityRadius(double impurityRadius) noexcept
    {
        _impurityRadius = asDerived()->checkImpurityRadius(impurityRadius);
    }

    void
    setTemperatures(const std::vector<double>& gridTemperatures);

    void
    syncClusterDataOnHost()
    {
        _subpaving.syncTiles(plsm::onHost);
        auto mirror = ClusterDataMirror(_subpaving, _gridSize);
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

    KOKKOS_INLINE_FUNCTION
    decltype(auto)
    getReactionRates(std::size_t reactionId)
    {
        return Kokkos::subview(_reactionRates, reactionId, Kokkos::ALL);
    }

    KOKKOS_INLINE_FUNCTION
    decltype(auto)
    getReactionCoefficients(std::size_t reactionId)
    {
        if (reactionId < _productionCoeffs.extent(0)) {
            return Kokkos::subview(_productionCoeffs, reactionId, Kokkos::ALL,
                Kokkos::ALL, Kokkos::ALL, Kokkos::ALL);
        }
        else {
            reactionId -= _productionCoeffs.extent(0);
            return Kokkos::subview(_dissociationCoeffs, reactionId, Kokkos::ALL,
                Kokkos::ALL, Kokkos::ALL, Kokkos::ALL);
        }
    }

    void
    computeAllFluxes(ConcentrationsView concentrations, FluxesView fluxes,
            std::size_t gridIndex);

    void
    computeAllPartials(ConcentrationsView concentrations,
            Kokkos::View<double*> values, std::size_t gridIndex);

    /**
     * Get the diagonal fill for the Jacobian, corresponding to the reactions.
     * Also populates the inverse map.
     *
     * @param fillMap Connectivity map.
     * @return The total number of partials.
     */
    std::size_t
    getDiagonalFill(SparseFillMap& fillMap);

protected:
    struct ClusterSetsPair
    {
        using ClusterSet = typename ReactionType::ClusterSet;
        std::vector<ClusterSet> prodClusterSets;
        std::vector<ClusterSet> dissClusterSets;
    };

    struct ClusterSetsViewPair
    {
        using ClusterSet = typename ReactionType::ClusterSet;
        Kokkos::View<ClusterSet*> prodClusterSets;
        Kokkos::View<ClusterSet*> dissClusterSets;
    };

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
        return _inverseMap(rowId, colId);
    }

private:
    Subpaving _subpaving;

    double _latticeParameter {};
    double _atomicVolume {};
    double _interstitialBias {};
    double _impurityRadius {};

    std::size_t _gridSize {};

    ClusterData _clusterData;
    ClusterDataMirror _clusterDataMirror;

    std::size_t _numDOFs {};

    Kokkos::View<ReactionType*> _reactions;
    Kokkos::View<double*****> _productionCoeffs;
    Kokkos::View<double*****> _dissociationCoeffs;
    Kokkos::View<double**> _reactionRates;

    // TODO: the original code uses an actual map here because it is sparse
    //       Look into Kokkos::Crs
    Kokkos::View<std::size_t**> _inverseMap;

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
    using ClusterData = typename Network::ClusterData;
    using ClusterDataRef = typename Network::ClusterDataRef;

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
    defineMomentIds(std::size_t& numDOFs);

    void
    defineReactions();

    std::size_t
    getDiagonalFill(typename Network::SparseFillMap& fillMap);
};
}









template <typename TReactionNetwork>
TReactionNetwork
makeReactionNetwork(
    const std::vector<typename TReactionNetwork::AmountType>& maxSpeciesAmounts,
    std::size_t gridSize, const IOptions& options)
{
    using Subpaving = typename TReactionNetwork::Subpaving;
    using Region = typename TReactionNetwork::Region;
    using Ival = typename TReactionNetwork::Ival;
    using AmountType = typename TReactionNetwork::AmountType;
    using SubdivRatio = typename Subpaving::SubdivisionRatioType;

    constexpr auto numSpecies = TReactionNetwork::getNumberOfSpecies();

    Region latticeRegion;
    SubdivRatio ratio;
    for (std::size_t i = 0; i < numSpecies; ++i) {
        auto maxAmt = maxSpeciesAmounts[i];
        latticeRegion[i] = Ival{0, maxAmt + 1};
        ratio[i] = maxAmt + 1;
    }
    Subpaving subpaving(latticeRegion, {ratio});

    //TODO: refine

    TReactionNetwork network(std::move(subpaving), gridSize, options);

    return network;
}


template <typename TReactionNetwork>
TReactionNetwork
makeSimpleReactionNetwork(
    typename TReactionNetwork::AmountType maxSpeciesAmount = 10)
{
    using Subpaving = typename TReactionNetwork::Subpaving;
    using Region = typename TReactionNetwork::Region;
    using Ival = typename TReactionNetwork::Ival;
    using AmountType = typename TReactionNetwork::AmountType;
    using SubdivRatio = typename Subpaving::SubdivisionRatioType;

    constexpr auto numSpecies = TReactionNetwork::getNumberOfSpecies();

    Ival ival{0, maxSpeciesAmount + 1};
    Region latticeRegion;
    SubdivRatio ratio;
    for (std::size_t i = 0; i < numSpecies; ++i) {
        latticeRegion[i] = ival;
        ratio[i] = maxSpeciesAmount + 1;
    }
    Subpaving subpaving(latticeRegion, {ratio});

    subpaving.refine(
        plsm::refine::RegionDetector<AmountType, numSpecies, plsm::Select>{
            latticeRegion});

    TReactionNetwork network(std::move(subpaving), 0);

    return network;
}
}
}

#include <experimental/ReactionNetwork.inl>
