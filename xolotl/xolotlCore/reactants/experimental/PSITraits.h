#pragma once

namespace xolotlCore
{
namespace experimental
{
template <typename TSpeciesEnum>
class PSIReactionNetwork;

template <typename TSpeciesEnum>
class PSIReaction;

template <typename TSpeciesEnum>
class PSIClusterGenerator;

enum class PSIFullSpeciesList
{
    He,
    D,
    T,
    V,
    I
};

template <>
struct HasInterstitial<PSIFullSpeciesList> : std::true_type
{
};

template <typename TSpeciesEnum>
struct ReactionNetworkTraits<PSIReactionNetwork<TSpeciesEnum>>
{
    using Species = TSpeciesEnum;

    static constexpr std::size_t numSpecies = 5;

    using ReactionType = PSIReaction<Species>;

    using ClusterGenerator = PSIClusterGenerator<Species>;
};
}
}
