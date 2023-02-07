#pragma once

#include <tuple>

#include <xolotl/core/network/ReactionNetworkTraits.h>

namespace xolotl
{
namespace core
{
namespace network
{
class NEProductionReaction;
class NEDissociationReaction;
class NEFullReSolutionReaction;
class NEPartialReSolutionReaction;
class NENucleationReaction;
class NEReactionNetwork;
class NEClusterGenerator;
namespace detail
{
class NEClusterUpdater;
}

enum class NESpecies
{
	Xe
};

inline const std::string&
toLabelString(NESpecies species)
{
	static const std::string labelArray[] = {"Xe"};
	return labelArray[static_cast<int>(species)];
}

inline const std::string&
toNameString(NESpecies species)
{
	static const std::string nameArray[] = {"Xenon"};
	return nameArray[static_cast<int>(species)];
}

template <>
struct NumberOfInterstitialSpecies<NESpecies> :
	std::integral_constant<std::size_t, 0>
{
};

template <>
struct NumberOfVacancySpecies<NESpecies> :
	std::integral_constant<std::size_t, 0>
{
};

template <>
struct ReactionNetworkTraits<NEReactionNetwork>
{
	using Species = NESpecies;

	static constexpr std::size_t numSpecies = 1;

	// using ReactionType = NEReaction;
	using ProductionReactionType = NEProductionReaction;
	using DissociationReactionType = NEDissociationReaction;
	using FullReSolutionReactionType = NEFullReSolutionReaction;
	using PartialReSolutionReactionType = NEPartialReSolutionReaction;
	using NucleationReactionType = NENucleationReaction;

	using ReactionTypeList = std::tuple<ProductionReactionType,
		DissociationReactionType, FullReSolutionReactionType,
		PartialReSolutionReactionType, NucleationReactionType>;

	using ClusterGenerator = NEClusterGenerator;
	using ClusterUpdater = detail::NEClusterUpdater;
};
} // namespace network
} // namespace core
} // namespace xolotl
