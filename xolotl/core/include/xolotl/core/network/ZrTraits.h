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
class ZrReactionNetwork;
class ZrProductionReaction;
class ZrDissociationReaction;
class ZrSinkReaction;
class ZrClusterGenerator;
namespace detail
{
class ZrClusterUpdater;
}

enum class ZrSpecies
{
	V,
	I
};

inline const std::string&
toLabelString(ZrSpecies species)
{
	static const std::string labelArray[] = {"V", "I"};
	return labelArray[static_cast<int>(species)];
}

inline const std::string&
toNameString(ZrSpecies species)
{
	static const std::string nameArray[] = {"Vacancy", "Interstitial"};
	return nameArray[static_cast<int>(species)];
}

template <>
struct NumberOfInterstitialSpecies<ZrSpecies> :
	std::integral_constant<std::size_t, 1>
{
};

template <>
struct NumberOfVacancySpecies<ZrSpecies> :
	std::integral_constant<std::size_t, 1>
{
};

template <>
struct SpeciesForGrouping<ZrSpecies, 2>
{
	using Sequence = EnumSequence<ZrSpecies, 2>;
	static constexpr auto first = Sequence(ZrSpecies::V);
	static constexpr auto last = Sequence(ZrSpecies::I);

	KOKKOS_INLINE_FUNCTION
	static constexpr std::underlying_type_t<ZrSpecies> mapToMomentId(
		EnumSequence<ZrSpecies, 2>)
	{
		return 0;
	}
};

template <>
struct ReactionNetworkTraits<ZrReactionNetwork>
{
	using Species = ZrSpecies;

	static constexpr std::size_t numSpecies = 2;

	using ProductionReactionType = ZrProductionReaction;
	using DissociationReactionType = ZrDissociationReaction;
	using SinkReactionType = ZrSinkReaction;

	using ReactionTypeList = std::tuple<ProductionReactionType,
		DissociationReactionType, SinkReactionType>;

	using ClusterGenerator = ZrClusterGenerator;
	using ClusterUpdater = detail::ZrClusterUpdater;
};
} // namespace network
} // namespace core
} // namespace xolotl
