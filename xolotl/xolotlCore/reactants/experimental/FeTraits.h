#pragma once

namespace xolotlCore
{
namespace experimental
{
template <typename TSpeciesEnum>
class FeReactionNetwork;

template <typename TSpeciesEnum>
class FeProductionReaction;

template <typename TSpeciesEnum>
class FeDissociationReaction;

template <typename TSpeciesEnum>
class FeClusterGenerator;

enum class FeFullSpeciesList
{
    He,
    V,
    I
};

template <>
struct HasInterstitial<FeFullSpeciesList> : std::true_type
{
};

template <typename TSpeciesEnum>
struct ReactionNetworkTraits<FeReactionNetwork<TSpeciesEnum>>
{
    using Species = TSpeciesEnum;

    static constexpr std::size_t numSpecies = 3;

    // using ReactionType = FeReaction<Species>;
    using ProductionReactionType = FeProductionReaction<Species>;
    using DissociationReactionType = FeDissociationReaction<Species>;

    using ClusterGenerator = FeClusterGenerator<Species>;
};
}
}
