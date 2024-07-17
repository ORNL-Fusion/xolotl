#pragma once

#include <xolotl/core/network/ReactionNetworkTraits.h>

namespace xolotl
{
namespace core
{
namespace network
{
class T91ReactionNetwork;
class T91ProductionReaction;
class T91DissociationReaction;
class T91SinkReaction;
class T91ClusterGenerator;

enum class T91SpeciesList
{
	He,
	V,
	I
};

inline const std::string&
toLabelString(T91SpeciesList species)
{
	static const std::string labelArray[] = {"He", "V", "I"};
	return labelArray[static_cast<int>(species)];
}

inline const std::string&
toNameString(T91SpeciesList species)
{
	static const std::string nameArray[] = {
		"Helium", "Vacancy", "Interstitial"};
	return nameArray[static_cast<int>(species)];
}

template <>
struct NumberOfInterstitialSpecies<T91SpeciesList> :
	std::integral_constant<std::size_t, 1>
{
};

template <>
struct NumberOfVacancySpecies<T91SpeciesList> :
	std::integral_constant<std::size_t, 1>
{
};

template <>
struct ReactionNetworkTraits<T91ReactionNetwork>
{
	using Species = T91SpeciesList;

	static constexpr std::size_t numSpecies = 3;

	using ProductionReactionType = T91ProductionReaction;
	using DissociationReactionType = T91DissociationReaction;
	using SinkReactionType = T91SinkReaction;

	using ReactionTypeList = std::tuple<ProductionReactionType,
		DissociationReactionType, SinkReactionType>;

	using ClusterGenerator = T91ClusterGenerator;
};
} // namespace network
} // namespace core
} // namespace xolotl
