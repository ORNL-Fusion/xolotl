#pragma once

namespace xolotlCore
{
namespace experimental
{
class FeReactionNetwork;
class FeProductionReaction;
class FeDissociationReaction;
class FeSinkReaction;
class FeClusterGenerator;

enum class FeSpeciesList
{
    He,
    V,
    I
};

template <>
struct HasInterstitial<FeSpeciesList> : std::true_type
{
};

template <>
struct ReactionNetworkTraits<FeReactionNetwork>
{
    using Species = FeSpeciesList;

    static constexpr std::size_t numSpecies = 3;

    using ProductionReactionType = FeProductionReaction;
    using DissociationReactionType = FeDissociationReaction;
    using SinkReactionType = FeSinkReaction;

    using ClusterGenerator = FeClusterGenerator;
};
}
}
