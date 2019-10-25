#pragma once

#include <cstddef>
#include <cstdint>
#include <type_traits>

#include <plsm/Subpaving.h>
#include <plsm/refine/RegionDetector.h>

#include <experimental/Reaction.h>
#include <experimental/SequencedEnumBase.h>

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

private:
    static constexpr std::size_t numSpecies = Traits::numSpecies;

public:
    using Species = typename Traits::Species;
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


class PSIReactionNetwork;


template <>
struct ReactionNetworkTraits<PSIReactionNetwork>
{
    enum class Species
    {
        V,
        I,
        He,
        D,
        T,
        first = V,
        last = T
    };


    // struct Species : SequencedEnumBase<Species, 5>
    // {
    //     using S = SequencedEnumBase<Species, 5>;
    //     using S::S;

    //     static constexpr S He{0};
    //     static constexpr S D{1};
    //     static constexpr S T{2};
    //     static constexpr S V{3};
    //     static constexpr S I{4};
    // };

    static constexpr std::size_t numSpecies = 5;

    using ReactionType = PSIReaction;
};


class PSIReactionNetwork : public ReactionNetwork<PSIReactionNetwork>
{
    using ReactionNetwork<PSIReactionNetwork>::ReactionNetwork;
};
}
}

#include <experimental/Cluster.h>
