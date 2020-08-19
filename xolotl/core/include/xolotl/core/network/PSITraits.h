#pragma once

#include <tuple>

#include <xolotl/core/network/ReactionNetworkTraits.h>

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

inline const std::string&
toLabelString(PSIFullSpeciesList species)
{
	static const std::string labelArray[] = {"He", "D", "T", "V", "I"};
	return labelArray[static_cast<int>(species)];
}

inline const std::string&
toNameString(PSIFullSpeciesList species)
{
	static const std::string nameArray[] = {
		"Helium", "Deuterium", "Tritium", "Vacancy", "Interstitial"};
	return nameArray[static_cast<int>(species)];
}

template <>
struct NumberOfInterstitialSpecies<PSIFullSpeciesList> :
	std::integral_constant<std::size_t, 1>
{
};

template <>
struct NumberOfVacancySpecies<PSIFullSpeciesList> :
	std::integral_constant<std::size_t, 1>
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
} // namespace network
} // namespace core
} // namespace xolotl
