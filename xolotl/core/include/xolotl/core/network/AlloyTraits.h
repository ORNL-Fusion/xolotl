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
class AlloyReactionNetwork;
class AlloyProductionReaction;
class AlloyDissociationReaction;
class AlloySinkReaction;
class AlloyTransformReaction;
class AlloyClusterGenerator;

enum class AlloySpecies
{
	V,
	PerfectV,
	FaultedV,
	I,
	PerfectI,
	FaultedI
};

inline const std::string&
toLabelString(AlloySpecies species)
{
	static const std::string labelArray[] = {
		"V", "PerfectV", "FaultedV", "I", "PerfectI", "FaultedI"};
	return labelArray[static_cast<int>(species)];
}

inline const std::string&
toNameString(AlloySpecies species)
{
	static const std::string nameArray[] = {"Vacancy", "PerfectV", "FaultedV",
		"Interstitial", "PerfectI", "FaultedI"};
	return nameArray[static_cast<int>(species)];
}

template <>
struct NumberOfInterstitialSpecies<AlloySpecies> :
	std::integral_constant<std::size_t, 3>
{
};

template <>
struct NumberOfVacancySpecies<AlloySpecies> :
	std::integral_constant<std::size_t, 3>
{
};

template <>
struct SpeciesForGrouping<AlloySpecies, 6>
{
	using Sequence = EnumSequence<AlloySpecies, 6>;
	static constexpr auto first = Sequence(AlloySpecies::V);
	static constexpr auto last = Sequence(AlloySpecies::FaultedI);

	KOKKOS_INLINE_FUNCTION
	static constexpr std::underlying_type_t<AlloySpecies>
	mapToMomentId(EnumSequence<AlloySpecies, 6>)
	{
		return 0;
	}
};

template <>
struct ReactionNetworkTraits<AlloyReactionNetwork>
{
	using Species = AlloySpecies;

	static constexpr std::size_t numSpecies = 6;

	using ProductionReactionType = AlloyProductionReaction;
	using DissociationReactionType = AlloyDissociationReaction;
	using SinkReactionType = AlloySinkReaction;
	using TransformReactionType = AlloyTransformReaction;

	using ReactionTypeList = std::tuple<ProductionReactionType,
		DissociationReactionType, SinkReactionType, TransformReactionType>;

	using ClusterGenerator = AlloyClusterGenerator;
};
} // namespace network
} // namespace core
} // namespace xolotl
