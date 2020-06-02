#pragma once

namespace xolotl
{
namespace core
{
namespace network
{
class AlloyReactionNetwork;
class AlloyProductionReaction;
class AlloyDissociationReaction;
class AlloySinkReaction;
// class AlloyReSolutionReaction;
class AlloyClusterGenerator;

enum class AlloySpecies
{
    V,
    Void,
    Faulted,
    I,
    Perfect,
    Frank
};

inline const char* toString(AlloySpecies specie) {
    static const char* nameArray[] = {"V", "Void", "Faulted", "I", "Perfect", "Frank"};
    return nameArray[static_cast<int>(specie)];
}

template <>
struct NumberOfInterstitialSpecies<AlloySpecies> : std::integral_constant<std::size_t,3>
{
};

template <>
struct NumberOfVacancySpecies<AlloySpecies> : std::integral_constant<std::size_t,3>
{
};

template <>
struct ReactionNetworkTraits<AlloyReactionNetwork>
{
    using Species = AlloySpecies;

    static constexpr std::size_t numSpecies = 6;

    using ProductionReactionType = AlloyProductionReaction;
    using DissociationReactionType = AlloyDissociationReaction;
    using SinkReactionType = AlloySinkReaction;

    using ReactionTypeList =
        std::tuple<ProductionReactionType, DissociationReactionType,
            SinkReactionType>;

    using ClusterGenerator = AlloyClusterGenerator;
};
}
}
}
