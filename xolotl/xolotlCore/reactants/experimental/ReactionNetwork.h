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
        return
            SpeciesRange(SpeciesSequence::first(), SpeciesSequence::lastNoI());
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

private:
    Subpaving _subpaving;
    Kokkos::View<ReactionType> _reactions;
};


template <typename TImpl>
ReactionNetwork<TImpl>::ReactionNetwork(AmountType maxSpeciesAmount)
    :
    _subpaving(Region{{
        Ival{0, maxSpeciesAmount+1},
        Ival{0, maxSpeciesAmount+1},
        Ival{0, maxSpeciesAmount+1},
        Ival{0, maxSpeciesAmount+1},
        Ival{0, maxSpeciesAmount+1}}},
    {{{
        maxSpeciesAmount+1,
        maxSpeciesAmount+1,
        maxSpeciesAmount+1,
        maxSpeciesAmount+1,
        maxSpeciesAmount+1}}})
{
}


template <typename TImpl>
typename ReactionNetwork<TImpl>::Cluster
ReactionNetwork<TImpl>::getCluster(const Composition& comp)
{
    _subpaving.syncAll(plsm::onHost);
    Cluster ret(*this, _subpaving.findTileId(comp, plsm::onHost));
    return ret;
}


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

    return network;
}


template <typename TSpeciesEnum>
class PSIReactionNetwork;


enum class PSIFullSpeciesList
{
    He,
    D,
    T,
    V,
    I
};


template <>
struct hasInterstitial<PSIFullSpeciesList> : std::true_type { };


template <typename TSpeciesEnum>
struct ReactionNetworkTraits<PSIReactionNetwork<TSpeciesEnum>>
{
    using Species = TSpeciesEnum;

    static constexpr std::size_t numSpecies = 5;

    using ReactionType = PSIReaction;
};


template <typename TSpeciesEnum>
class PSIReactionNetwork :
    public ReactionNetwork<PSIReactionNetwork<TSpeciesEnum>>
{
public:
    using ReactionNetwork<PSIReactionNetwork<TSpeciesEnum>>::ReactionNetwork;
};
}
}

#include <experimental/Cluster.h>
