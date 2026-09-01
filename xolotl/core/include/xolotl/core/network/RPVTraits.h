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
	S,
	V,
	I,
	Loop
};

inline const std::string&
toLabelString(RPVSpeciesList species)
{
	static const std::string labelArray[] = {"S", "V", "I", "Loop"};
	return labelArray[static_cast<int>(species)];
}

inline const std::string&
toNameString(RPVSpeciesList species)
{
	static const std::string nameArray[] = {
		"Solute", "Vacancy", "Interstitial", "Loop"};
	return nameArray[static_cast<int>(species)];
}

template <>
struct NumberOfInterstitialSpecies<RPVSpeciesList> :
	std::integral_constant<std::size_t, 2>
{
};

template <>
struct NumberOfVacancySpecies<RPVSpeciesList> :
	std::integral_constant<std::size_t, 1>
{
};

template <>
struct SpeciesForGrouping<RPVSpeciesList, 4>
{
	using Sequence = EnumSequence<RPVSpeciesList, 4>;
	static constexpr auto first = Sequence(RPVSpeciesList::S);
	static constexpr auto last = Sequence(RPVSpeciesList::Loop);

	KOKKOS_INLINE_FUNCTION
	static constexpr std::underlying_type_t<RPVSpeciesList>
	mapToMomentId(EnumSequence<RPVSpeciesList, 4> value)
	{
		if (value == RPVSpeciesList::I or value == RPVSpeciesList::Loop)
			return 1;
		return value();
	}
};

template <>
struct ReactionNetworkTraits<RPVReactionNetwork>
{
	using Species = RPVSpeciesList;

	static constexpr std::size_t numSpecies = 4;

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
