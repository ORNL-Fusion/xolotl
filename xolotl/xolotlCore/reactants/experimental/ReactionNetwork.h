#pragma once

#include <cstddef>
#include <cstdint>
#include <type_traits>

#include <plsm/Subpaving.h>
#include <plsm/refine/RegionDetector.h>

#include <Species.h>

namespace xolotlCore
{
namespace experimental
{
template <std::size_t NumSpecies>
class ReactionNetwork
{
public:
    //TODO: Maybe we could traits-ify this:
    //Define a phase space for a specific application.
    //The "phase space" could define the underlying amount type to
    //use and the "Species" enum and/or number of species

    using AmountType = std::uint32_t;
    using Subpaving = plsm::Subpaving<AmountType, NumSpecies, Species>;
    using Composition = typename Subpaving::PointType;
    using Region = typename Subpaving::RegionType;
    using Ival = typename Region::IntervalType;

    class Cluster;

    ReactionNetwork() = delete;
    ReactionNetwork(AmountType maxSpeciesAmount);

    Cluster
    get(const Composition& comp);

    Subpaving&
    getSubpaving()
    {
        return _subpaving;
    }

private:
    Subpaving _subpaving;
};


template <std::size_t NumSpecies>
ReactionNetwork<NumSpecies>::ReactionNetwork(AmountType maxSpeciesAmount)
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


template <std::size_t NumSpecies>
typename ReactionNetwork<NumSpecies>::Cluster
ReactionNetwork<NumSpecies>::get(const Composition& comp)
{
    Cluster ret(*this, _subpaving.getTileId(comp));
    return ret;
}


template <std::size_t NumSpecies>
ReactionNetwork<NumSpecies>
makeSimpleReactionNetwork(
    typename ReactionNetwork<NumSpecies>::AmountType maxSpeciesAmount = 10)
{
    using AmountType = typename ReactionNetwork<NumSpecies>::AmountType;
    ReactionNetwork<NumSpecies> network(maxSpeciesAmount);

    network.getSubpaving().refine(
        plsm::refine::RegionDetector<AmountType, NumSpecies, plsm::Select>{
            network.getSubpaving().getLatticeRegion()});

    return network;
}
}
}

#include <experimental/Cluster.h>
