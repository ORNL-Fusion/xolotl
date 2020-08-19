#pragma once

#include <xolotl/core/network/ReactionNetworkTraits.h>

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
class FeClusterGenerator;

enum class FeSpeciesList
{
	He,
	V,
	I
};

inline const std::string&
toLabelString(FeSpeciesList species)
{
	static const std::string labelArray[] = {"He", "V", "I"};
	return labelArray[static_cast<int>(species)];
}

inline const std::string&
toNameString(FeSpeciesList species)
{
	static const std::string nameArray[] = {
		"Helium", "Vacancy", "Interstitial"};
	return nameArray[static_cast<int>(species)];
}

template <>
struct NumberOfInterstitialSpecies<FeSpeciesList> :
	std::integral_constant<std::size_t, 1>
{
};

template <>
struct NumberOfVacancySpecies<FeSpeciesList> :
	std::integral_constant<std::size_t, 1>
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

	using ReactionTypeList = std::tuple<ProductionReactionType,
		DissociationReactionType, SinkReactionType>;

	using ClusterGenerator = FeClusterGenerator;
};
} // namespace network
} // namespace core
} // namespace xolotl
