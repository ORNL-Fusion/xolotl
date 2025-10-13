#pragma once

#include <xolotl/core/network/ReactionNetworkTraits.h>

namespace xolotl
{
namespace core
{
namespace network
{
class RPVReactionNetwork;
class RPVProductionReaction;
class RPVDissociationReaction;
class RPVSinkReaction;
class RPVClusterGenerator;

enum class RPVSpeciesList
{
	V,
	I
};

inline const std::string&
toLabelString(RPVSpeciesList species)
{
	static const std::string labelArray[] = {"V", "I"};
	return labelArray[static_cast<int>(species)];
}

inline const std::string&
toNameString(RPVSpeciesList species)
{
	static const std::string nameArray[] = {"Vacancy", "Interstitial"};
	return nameArray[static_cast<int>(species)];
}

template <>
struct NumberOfInterstitialSpecies<RPVSpeciesList> :
	std::integral_constant<std::size_t, 1>
{
};

template <>
struct NumberOfVacancySpecies<RPVSpeciesList> :
	std::integral_constant<std::size_t, 1>
{
};

template <>
struct ReactionNetworkTraits<RPVReactionNetwork>
{
	using Species = RPVSpeciesList;

	static constexpr std::size_t numSpecies = 2;

	using ProductionReactionType = RPVProductionReaction;
	using DissociationReactionType = RPVDissociationReaction;
	using SinkReactionType = RPVSinkReaction;

	using ReactionTypeList = std::tuple<ProductionReactionType,
		DissociationReactionType, SinkReactionType>;

	using ClusterGenerator = RPVClusterGenerator;
};
} // namespace network
} // namespace core
} // namespace xolotl
