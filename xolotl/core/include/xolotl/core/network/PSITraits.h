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
template <typename TSpeciesEnum>
class PSIReactionNetwork;

template <typename TSpeciesEnum>
class PSIProductionReaction;

template <typename TSpeciesEnum>
class PSIDissociationReaction;

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
struct NumberOfInterstitialSpecies<PSIFullSpeciesList> : std::integral_constant<std::size_t,1>
{
};

template <>
struct NumberOfVacancySpecies<PSIFullSpeciesList> : std::integral_constant<std::size_t,1>
{
};

template <typename TSpeciesEnum>
struct ReactionNetworkTraits<PSIReactionNetwork<TSpeciesEnum>>
{
    using Species = TSpeciesEnum;

    static constexpr std::size_t numSpecies = 5;

    using ProductionReactionType = PSIProductionReaction<Species>;
    using DissociationReactionType = PSIDissociationReaction<Species>;

    using ReactionTypeList =
        std::tuple<ProductionReactionType, DissociationReactionType>;

    using ClusterGenerator = PSIClusterGenerator<Species>;
};
}
}
}
