#pragma once

namespace xolotl
{
namespace core
{
namespace network
{
class FeReactionNetwork;
class FeProductionReaction;
class FeDissociationReaction;
class FeSinkReaction;
// class FeReSolutionReaction;
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
struct NumberOfInterstitialSpecies<FeSpeciesList> : std::integral_constant<std::size_t,1>
{
};

template <>
struct NumberOfVacancySpecies<FeSpeciesList> : std::integral_constant<std::size_t,1>
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

    using ReactionTypeList =
        std::tuple<ProductionReactionType, DissociationReactionType,
            SinkReactionType>;

    using ClusterGenerator = FeClusterGenerator;
};
}
}
}
