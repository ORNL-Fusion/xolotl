#pragma once

#include <cstddef>
#include <cstdint>
#include <type_traits>

#include <plsm/Subpaving.h>
#include <plsm/refine/RegionDetector.h>

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
public:
    using Traits = ReactionNetworkTraits<TImpl>;
    using Species = typename Traits::Species;

private:
    static constexpr std::size_t numSpecies = Traits::numSpecies;

public:
    using SpeciesSequence = SpeciesEnumSequence<Species, numSpecies>;
    using SpeciesRange = EnumSequenceRange<Species, numSpecies>;
    using ReactionType = typename Traits::ReactionType;
    using AmountType = std::uint32_t;
    using Subpaving = plsm::Subpaving<AmountType, numSpecies, Species>;
    using Composition = typename Subpaving::PointType;
    using Region = typename Subpaving::RegionType;
    using Ival = typename Region::IntervalType;

    class Cluster;

    template <typename T>
    class Reaction;

    ReactionNetwork() = delete;

    //TODO: Need a more versatile constructor interface
    //      (and probably don't need the 'make' function)
    ReactionNetwork(AmountType maxSpeciesAmount);

    static
    constexpr std::size_t
    getNumberOfSpecies() noexcept
    {
        return numSpecies;
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
    setLatticeParameter(double latticeParameter)
    {
        _latticeParameter = latticeParameter;
        _atomicVolume =
            0.5 * latticeParameter * latticeParameter * latticeParameter;
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
    getMomentIds(std::size_t clusterId)
    {
        return Kokkos::subview(_momentIds, clusterId, Kokkos::ALL);
    }

    double
    getReactionRadius(std::size_t clusterId)
    {
        return _reactionRadius(clusterId);
    }

    void
    defineMomentIds();

    decltype(auto)
    getReactionRates(std::size_t reactionId)
    {
        return Kokkos::subview(_reactionRates, reactionId, Kokkos::ALL);
    }

private:
    double
    getTemperature(std::size_t gridIndex) const noexcept
    {
        return _temperature(gridIndex);
    }

    double
    getDiffusionCoefficient(std::size_t reactionId, std::size_t gridIndex)
    {
        return _diffusionCoefficient(reactionId, gridIndex);
    }

private:
    double _latticeParameter;
    double _atomicVolume;
    Kokkos::View<double*> _temperature;

    Subpaving _subpaving;

    Kokkos::View<std::size_t*[4]> _momentIds;
    Kokkos::View<double*> _reactionRadius;
    Kokkos::View<double**> _diffusionCoefficient;
    Kokkos::View<double*> _formationEnergy;

    Kokkos::View<ReactionType*> _reactions;
    Kokkos::View<double**> _reactionRates;
};


template <typename TReactionNetwork>
TReactionNetwork
makeSimpleReactionNetwork(
    typename TReactionNetwork::AmountType maxSpeciesAmount = 10)
{
    using AmountType = typename TReactionNetwork::AmountType;
    TReactionNetwork network(maxSpeciesAmount);

    constexpr auto numSpecies = network.getNumberOfSpecies();
    network.getSubpaving().refine(
        plsm::refine::RegionDetector<AmountType, numSpecies, plsm::Select>{
            network.getSubpaving().getLatticeRegion()});

    network.defineMomentIds();

    return network;
}
}
}

#include <experimental/ReactionNetwork.inl>
#include <experimental/Cluster.h>
#include <experimental/Reaction.h>
