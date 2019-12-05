#pragma once

#include <cstddef>
#include <cstdint>
#include <type_traits>

#include <plsm/Subpaving.h>
#include <plsm/refine/RegionDetector.h>

#include <IOptions.h>
#include <Options.h>

#include <experimental/EnumSequence.h>

namespace xolotlCore
{
namespace experimental
{
template <typename TImpl>
struct ReactionNetworkTraits
{
};


template <typename TImpl>
class ReactionNetwork
{
    static constexpr auto invalid = plsm::invalid<std::size_t>;

public:
    using Traits = ReactionNetworkTraits<TImpl>;
    using Species = typename Traits::Species;

private:
    static constexpr std::size_t numSpecies = Traits::numSpecies;

public:
    using SpeciesSequence = SpeciesEnumSequence<Species, numSpecies>;
    using SpeciesRange = EnumSequenceRange<Species, numSpecies>;
    using ReactionType = typename Traits::ReactionType;
    using ClusterGenerator = typename Traits::ClusterGenerator;
    using AmountType = std::uint32_t;
    using Subpaving = plsm::Subpaving<AmountType, numSpecies, Species>;
    using SubdivisionRatio = typename Subpaving::SubdivisionRatioType;
    using Composition = typename Subpaving::PointType;
    using Region = typename Subpaving::RegionType;
    using Ival = typename Region::IntervalType;
    using ConcentrationsView = Kokkos::View<double*, Kokkos::MemoryUnmanaged>;
    using FluxesView = Kokkos::View<double*, Kokkos::MemoryUnmanaged>;
    using ConnectivityView = Kokkos::View<size_t**, Kokkos::MemoryUnmanaged>;
    using SparseFillMap = std::unordered_map<int, std::vector<int>>;

    class Cluster;

    template <typename T>
    class Reaction;

    ReactionNetwork() = delete;

    ReactionNetwork(Subpaving&& subpaving, std::size_t gridSize,
        const IOptions& options);

    ReactionNetwork(Subpaving&& subpaving, std::size_t gridSize);

    //TODO: Need a more versatile constructor interface
    //      (and probably don't need the 'make' function)

    ReactionNetwork(const std::vector<AmountType>& maxSpeciesAmounts,
        const std::vector<SubdivisionRatio>& subdivisionRatios,
        std::size_t gridSize, const IOptions& options);

    ReactionNetwork(const std::vector<AmountType>& maxSpeciesAmounts,
        std::size_t gridSize, const IOptions& options);

    static
    constexpr std::size_t
    getNumberOfSpecies() noexcept
    {
        return SpeciesSequence::size();
    }

    static
    constexpr std::size_t
    getNumberOfSpeciesNoI() noexcept
    {
        return SpeciesSequence::sizeNoI();
    }

    static
    constexpr SpeciesRange
    getSpeciesRange() noexcept
    {
        return SpeciesRange{};
    }

    static
    constexpr SpeciesRange
    getSpeciesRangeNoI() noexcept
    {
        return SpeciesRange(SpeciesSequence::first(),
            SpeciesSequence::lastNoI());
    }

    std::size_t
    getDOF() const noexcept
    {
        return _numDOFs;
    }

    double
    getLatticeParameter() const noexcept
    {
        return _latticeParameter;
    }

    double
    getAtomicVolume() const noexcept
    {
        return _atomicVolume;
    }

    void
    setLatticeParameter(double latticeParameter);

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

    Cluster
    findCluster(const Composition& comp)
    {
        //FIXME: explicitly using host space
        _subpaving.syncAll(plsm::onHost);
        return Cluster(*this, _subpaving.findTileId(comp, plsm::onHost));
    }

    Cluster
    getCluster(std::size_t clusterId)
    {
        return Cluster(*this, clusterId);
    }

    Subpaving&
    getSubpaving()
    {
        return _subpaving;
    }

    decltype(auto)
    getReactionRates(std::size_t reactionId)
    {
        return Kokkos::subview(_reactionRates, reactionId, Kokkos::ALL);
    }

    void
    computeAllFluxes(ConcentrationsView concentrations, FluxesView fluxes,
            std::size_t gridIndex);

    void
    computeAllPartials(ConcentrationsView concentrations,
            Kokkos::View<double*> values, std::size_t gridIndex);

private:
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

    double
    getTemperature(std::size_t gridIndex) const noexcept
    {
        return _temperature(gridIndex);
    }

    size_t
    getReactionIndex(std::size_t rowId, std::size_t colId) const noexcept
    {
        return _inverseMap(rowId, colId);
    }

    /**
     * Get the diagonal fill for the Jacobian, corresponding to the reactions.
     * Also populates the inverse map.
     *
     * @param fillMap Connectivity map.
     * @return The total number of partials.
     */
    size_t
    getDiagonalFill(SparseFillMap& fillMap);

private:
    Subpaving _subpaving;

    double _latticeParameter {};
    double _atomicVolume {};
    double _interstitialBias {};
    double _impurityRadius {};

    Kokkos::View<double*> _temperature;

    std::size_t _numClusters {};
    std::size_t _numDOFs {};
    Kokkos::View<std::size_t*[4]> _momentIds;
    Kokkos::View<double*> _reactionRadius;
    Kokkos::View<double**> _diffusionCoefficient;
    Kokkos::View<double*> _formationEnergy;
    Kokkos::View<double*> _migrationEnergy;
    Kokkos::View<double*> _diffusionFactor;

    Kokkos::View<ReactionType*> _reactions;

    Kokkos::View<double**> _reactionRates;

    // TODO: the original code uses an actual map here because it is sparse
    Kokkos::View<size_t**> _inverseMap;
};


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

#include <experimental/Cluster.h>
#include <experimental/Reaction.h>
#include <experimental/ReactionNetwork.inl>
