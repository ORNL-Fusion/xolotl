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
class FeCrReactionNetwork;
class FeCrProductionReaction;
class FeCrDissociationReaction;
class FeCrSinkReaction;
class FeCrTransformReaction;
class FeCrClusterGenerator;

enum class FeCrSpecies
{
	Trap,
	V,
	I,
	Free,
	Trapped,
	Junction,
	Complex,
	Loop
};

inline const std::string&
toLabelString(FeCrSpecies species)
{
	static const std::string labelArray[] = {
		"Trap", "V", "I", "Free", "Trapped", "Junction", "Complex", "Loop"};
	return labelArray[static_cast<int>(species)];
}

inline const std::string&
toNameString(FeCrSpecies species)
{
	static const std::string nameArray[] = {"Trap", "Void", "Interstitial110",
		"Free111", "Trapped111", "Junction", "Complex", "Loop100"};
	return nameArray[static_cast<int>(species)];
}

template <>
struct NumberOfInterstitialSpecies<FeCrSpecies> :
	std::integral_constant<std::size_t, 6>
{
};

template <>
struct NumberOfVacancySpecies<FeCrSpecies> :
	std::integral_constant<std::size_t, 1>
{
};

template <>
struct SpeciesForGrouping<FeCrSpecies, 8>
{
	using Sequence = EnumSequence<FeCrSpecies, 8>;
	static constexpr auto first = Sequence(FeCrSpecies::Trap);
	static constexpr auto last = Sequence(FeCrSpecies::Loop);

	KOKKOS_INLINE_FUNCTION
	static constexpr std::underlying_type_t<FeCrSpecies>
	mapToMomentId(EnumSequence<FeCrSpecies, 8>)
	{
		return 0;
	}
};

template <>
struct ReactionNetworkTraits<FeCrReactionNetwork>
{
	using Species = FeCrSpecies;

	static constexpr std::size_t numSpecies = 8;

	using ProductionReactionType = FeCrProductionReaction;
	using DissociationReactionType = FeCrDissociationReaction;
	using SinkReactionType = FeCrSinkReaction;
	using TransformReactionType = FeCrTransformReaction;

	using ReactionTypeList = std::tuple<ProductionReactionType,
		DissociationReactionType, SinkReactionType, TransformReactionType>;

	using ClusterGenerator = FeCrClusterGenerator;
};
} // namespace network
} // namespace core
} // namespace xolotl
