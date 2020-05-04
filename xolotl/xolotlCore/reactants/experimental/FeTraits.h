#pragma once

namespace xolotlCore
{
namespace experimental
{
class FeReactionNetwork;
class FeProductionReaction;
class FeDissociationReaction;
class FeSinkReaction;
class FeReSolutionReaction;
class FeClusterGenerator;

enum class FeSpeciesList
{
    He,
    V,
    I
};

inline const char* toString(FeSpeciesList specie) {
    static const char* nameArray[] = {"He", "V", "I"};
    return nameArray[static_cast<int>(specie)];
}

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
    using ReSolutionReactionType = FeReSolutionReaction;

    using ReactionTypeList =
        std::tuple<ProductionReactionType, DissociationReactionType,
            SinkReactionType, ReSolutionReactionType>;

    using ClusterGenerator = FeClusterGenerator;
};
}
}
