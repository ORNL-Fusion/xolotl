#pragma once

namespace xolotlCore
{
namespace experimental
{
class NEReaction;
class NEReactionNetwork;
class NEClusterGenerator;

enum class NESpecies
{
    Xe
};

template <>
struct HasInterstitial<NESpecies> : std::false_type
{
};

template <>
struct ReactionNetworkTraits<NEReactionNetwork>
{
    using Species = NESpecies;

    static constexpr std::size_t numSpecies = 1;

    using ReactionType = NEReaction;

    using ClusterGenerator = NEClusterGenerator;
};
}
}
