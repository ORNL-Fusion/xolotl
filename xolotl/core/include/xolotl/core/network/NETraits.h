#pragma once

#include <tuple>

#include <xolotl/core/network/ReactionNetworkTraits.h>
#include <xolotl/core/network/SpeciesEnumSequence.h>

namespace xolotl
{
namespace core
{
namespace network
{
class NEProductionReaction;
class NEDissociationReaction;
class NESinkReaction;
class NEReSolutionReaction;
class NEReactionNetwork;
class NEClusterGenerator;
namespace detail
{
class NEClusterUpdater;
}

enum class NESpecies
{
    Xe
};

inline const char* toString(NESpecies specie) {
    static const char* nameArray[] = {"Xe"};
    return nameArray[static_cast<int>(specie)];
}

template <>
struct NumberOfInterstitialSpecies<NESpecies> : std::integral_constant<std::size_t,0>
{
};

template <>
struct NumberOfVacancySpecies<NESpecies> : std::integral_constant<std::size_t,0>
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

    using ReactionTypeList =
        std::tuple<ProductionReactionType, DissociationReactionType,
            ReSolutionReactionType>;

    using ClusterGenerator = NEClusterGenerator;
    using ClusterUpdater = detail::NEClusterUpdater;
};
}
}
}
