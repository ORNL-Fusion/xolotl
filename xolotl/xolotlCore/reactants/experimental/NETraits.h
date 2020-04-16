#pragma once

namespace xolotlCore
{
namespace experimental
{
class NEProductionReaction;
class NEDissociationReaction;
class NESinkReaction;
class NEReSolutionReaction;
class NEReactionNetwork;
class NEClusterGenerator;

enum class NESpecies
{
    Xe
};

inline const char* toString(NESpecies specie) {
    static const char* nameArray[] = {"Xe"};
    return nameArray[static_cast<int>(specie)];
}

template <>
struct HasInterstitial<NESpecies> : std::false_type
{
};

template <>
struct ReactionNetworkTraits<NEReactionNetwork>
{
    using Species = NESpecies;

    static constexpr std::size_t numSpecies = 1;

    // using ReactionType = NEReaction;
    using ProductionReactionType = NEProductionReaction;
    using DissociationReactionType = NEDissociationReaction;
    using SinkReactionType = NESinkReaction;
    using ReSolutionReactionType = NEReSolutionReaction;

    using ClusterGenerator = NEClusterGenerator;
};
}
}
