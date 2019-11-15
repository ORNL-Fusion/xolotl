#pragma once

#include <experimental/ReactionNetwork.h>

namespace xolotlCore
{
namespace experimental
{
template <typename TSpeciesEnum>
class PSIReactionNetwork;


template <typename TSpeciesEnum>
class PSIReaction;


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

    using ReactionType = PSIReaction<TSpeciesEnum>;
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

#include <experimental/PSIReaction.h>
