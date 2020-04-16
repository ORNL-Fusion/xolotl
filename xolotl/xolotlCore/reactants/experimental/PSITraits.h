#pragma once

namespace xolotlCore
{
namespace experimental
{
template <typename TSpeciesEnum>
class PSIReactionNetwork;

template <typename TSpeciesEnum>
class PSIProductionReaction;

template <typename TSpeciesEnum>
class PSIDissociationReaction;

template <typename TSpeciesEnum>
class PSISinkReaction;

template <typename TSpeciesEnum>
class PSIReSolutionReaction;

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

inline const char* toString(PSIFullSpeciesList specie) {
    static const char* nameArray[] = {"He", "D", "T", "V", "I"};
    return nameArray[static_cast<int>(specie)];
}

template <>
struct HasInterstitial<PSIFullSpeciesList> : std::true_type
{
};

template <typename TSpeciesEnum>
struct ReactionNetworkTraits<PSIReactionNetwork<TSpeciesEnum>>
{
    using Species = TSpeciesEnum;

    static constexpr std::size_t numSpecies = 5;

    // using ReactionType = PSIReaction<Species>;
    using ProductionReactionType = PSIProductionReaction<Species>;
    using DissociationReactionType = PSIDissociationReaction<Species>;
    using SinkReactionType = PSISinkReaction<Species>;
    using ReSolutionReactionType = PSIReSolutionReaction<Species>;

    using ClusterGenerator = PSIClusterGenerator<Species>;
};
}
}
