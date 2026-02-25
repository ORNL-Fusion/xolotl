#pragma once

#include <xolotl/core/network/ReactionNetworkTraits.h>

namespace xolotl
{
namespace core
{
namespace network
{
class VReactionNetwork;
class VProductionReaction;
class VDissociationReaction;
class VSinkReaction;
class VClusterGenerator;

enum class VSpeciesList
{
	He,
	V,
	I
};

inline const std::string&
toLabelString(VSpeciesList species)
{
	static const std::string labelArray[] = {"He", "V", "I"};
	return labelArray[static_cast<int>(species)];
}

inline const std::string&
toNameString(VSpeciesList species)
{
	static const std::string nameArray[] = {
		"Helium", "Vacancy", "Interstitial"};
	return nameArray[static_cast<int>(species)];
}

template <>
struct NumberOfSpecies<VSpeciesList> : std::integral_constant<std::size_t, 3>
{
};

template <>
struct NumberOfInterstitialSpecies<VSpeciesList> :
	std::integral_constant<std::size_t, 1>
{
};

template <>
struct NumberOfVacancySpecies<VSpeciesList> :
	std::integral_constant<std::size_t, 1>
{
};

template <>
struct SpeciesForGrouping<VSpeciesList, 3>
{
	using Sequence = EnumSequence<VSpeciesList, 3>;
	static constexpr auto first = Sequence(VSpeciesList::He);
	static constexpr auto last = Sequence(VSpeciesList::I);

	KOKKOS_INLINE_FUNCTION
	static constexpr std::underlying_type_t<VSpeciesList>
	mapToMomentId(EnumSequence<VSpeciesList, 3> value)
	{
		if (value == VSpeciesList::I)
			return 1;
		return value();
	}
};

template <>
struct ReactionNetworkTraits<VReactionNetwork>
{
	using Species = VSpeciesList;

	static constexpr std::size_t numSpecies = 3;

	using ProductionReactionType = VProductionReaction;
	using DissociationReactionType = VDissociationReaction;
	using SinkReactionType = VSinkReaction;

	using ReactionTypeList = std::tuple<ProductionReactionType,
		DissociationReactionType, SinkReactionType>;

	using ClusterGenerator = VClusterGenerator;
};
} // namespace network
} // namespace core
} // namespace xolotl
