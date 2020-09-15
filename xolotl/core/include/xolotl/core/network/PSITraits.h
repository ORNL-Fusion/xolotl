#pragma once

#include <tuple>

#include <xolotl/core/network/ReactionNetworkTraits.h>

namespace xolotl
{
namespace core
{
namespace network
{
namespace detail
{
template <typename TSpeciesEnum>
struct HasDeuterium : std::false_type
{
};

template <typename TSpeciesEnum>
struct HasTritium : std::false_type
{
};
}

template <typename TSpeciesEnum>
inline constexpr bool hasDeuterium = detail::HasDeuterium<TSpeciesEnum>::value;

template <typename TSpeciesEnum>
inline constexpr bool hasTritium = detail::HasTritium<TSpeciesEnum>::value;

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

namespace detail
{
template <>
struct HasDeuterium<PSIFullSpeciesList> : std::true_type
{
};

template <>
struct HasTritium<PSIFullSpeciesList> : std::true_type
{
};
}

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
struct NumberOfSpecies<PSIFullSpeciesList> :
	std::integral_constant<std::size_t, 5>
{
};

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

enum class PSIHeliumSpeciesList
{
	He,
	V,
	I
};

inline const std::string&
toLabelString(PSIHeliumSpeciesList species)
{
	static const std::string labelArray[] = {"He", "V", "I"};
	return labelArray[static_cast<int>(species)];
}

inline const std::string&
toNameString(PSIHeliumSpeciesList species)
{
	static const std::string nameArray[] = {
		"Helium", "Vacancy", "Interstitial"};
	return nameArray[static_cast<int>(species)];
}

template <>
struct NumberOfSpecies<PSIHeliumSpeciesList> :
	std::integral_constant<std::size_t, 3>
{
};

template <>
struct NumberOfInterstitialSpecies<PSIHeliumSpeciesList> :
	std::integral_constant<std::size_t, 1>
{
};

template <>
struct NumberOfVacancySpecies<PSIHeliumSpeciesList> :
	std::integral_constant<std::size_t, 1>
{
};

template <typename TSpeciesEnum>
struct ReactionNetworkTraits<PSIReactionNetwork<TSpeciesEnum>>
{
	using Species = TSpeciesEnum;

	static constexpr std::size_t numSpecies = numberOfSpecies<TSpeciesEnum>();

	using ProductionReactionType = PSIProductionReaction<Species>;
	using DissociationReactionType = PSIDissociationReaction<Species>;

	using ReactionTypeList =
		std::tuple<ProductionReactionType, DissociationReactionType>;

	using ClusterGenerator = PSIClusterGenerator<Species>;
};
} // namespace network
} // namespace core
} // namespace xolotl
