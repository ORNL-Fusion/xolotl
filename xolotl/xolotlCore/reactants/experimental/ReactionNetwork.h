#pragma once

#include <cstddef>
#include <cstdint>
#include <type_traits>

#include <plsm/Subpaving.h>
#include <plsm/refine/RegionDetector.h>

#include <experimental/Reaction.h>
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

    Cluster
    getCluster(const Composition& comp);

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

    void
    defineMomentIds();

private:
    Subpaving _subpaving;

    Kokkos::View<std::size_t*[4]> _momentIds;
    Kokkos::View<double*> _reactionRadius;
    Kokkos::View<double**> _diffusionCoefficient;

    Kokkos::View<ReactionType*> _reactions;
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
