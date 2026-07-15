#pragma once

#include <tuple>

#include <Kokkos_UnorderedMap.hpp>

#include <xolotl/core/network/ReactionNetworkTraits.h>
#include <xolotl/core/network/detail/ClusterData.h>

namespace xolotl
{
namespace core
{
namespace network
{
class LiProductionReaction;
class LiDissociationReaction;
class LiTrapReaction;
class LiReactionNetwork;
class LiClusterGenerator;

enum class LiSpecies
{
        Trap,
	H,
	V
};

inline const std::string&
toLabelString(LiSpecies species)
{
	static const std::string labelArray[] = {"Trap", "H", "V"};
	return labelArray[static_cast<int>(species)];
}

inline const std::string&
toNameString(LiSpecies species)
{
	static const std::string nameArray[] = {"Trap", "Hydrogen", "Vacancy"};
	return nameArray[static_cast<int>(species)];
}

template <>
struct NumberOfSpecies<LiSpecies> : std::integral_constant<std::size_t, 3>
{
};

template <>
struct NumberOfInterstitialSpecies<LiSpecies> :
	std::integral_constant<std::size_t, 0>
{
};

template <>
struct NumberOfVacancySpecies<LiSpecies> :
	std::integral_constant<std::size_t, 1>
{
};

template <>
struct SpeciesForGrouping<LiSpecies, 3>
{
	using Sequence = EnumSequence<LiSpecies, 3>;
	static constexpr auto first = Sequence(LiSpecies::H);
	static constexpr auto last = Sequence(LiSpecies::V);

	KOKKOS_INLINE_FUNCTION
	static constexpr std::underlying_type_t<LiSpecies>
	mapToMomentId(EnumSequence<LiSpecies, 3> value)
	{
		return value() - 1;
	}
};

template <>
struct ReactionNetworkTraits<LiReactionNetwork>
{
	using Species = LiSpecies;

	static constexpr std::size_t numSpecies = 3;

	// using ReactionType = LiReaction;
	using ProductionReactionType = LiProductionReaction;
	using DissociationReactionType = LiDissociationReaction;
	using TrapReactionType = LiTrapReaction;

	using ReactionTypeList =
		std::tuple<ProductionReactionType, DissociationReactionType, TrapReactionType>;

	using ClusterGenerator = LiClusterGenerator;
};

} // namespace network
} // namespace core
} // namespace xolotl
