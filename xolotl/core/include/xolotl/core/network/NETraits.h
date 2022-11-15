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
class NEProductionReaction;
class NEDissociationReaction;
class NEReSolutionReaction;
class NENucleationReaction;
class NESinkReaction;
class NEReactionNetwork;
class NEClusterGenerator;
namespace detail
{
class NEClusterUpdater;
}

enum class NESpecies
{
	Xe,
	V,
	I
};

inline const std::string&
toLabelString(NESpecies species)
{
	static const std::string labelArray[] = {"Xe", "V", "I"};
	return labelArray[static_cast<int>(species)];
}

inline const std::string&
toNameString(NESpecies species)
{
	static const std::string nameArray[] = {"Xenon", "Vacancy", "Interstitial"};
	return nameArray[static_cast<int>(species)];
}

template <>
struct NumberOfInterstitialSpecies<NESpecies> :
	std::integral_constant<std::size_t, 1>
{
};

template <>
struct NumberOfVacancySpecies<NESpecies> :
	std::integral_constant<std::size_t, 1>
{
};

template <>
struct ReactionNetworkTraits<NEReactionNetwork>
{
	using Species = NESpecies;

	static constexpr std::size_t numSpecies = 3;

	// using ReactionType = NEReaction;
	using ProductionReactionType = NEProductionReaction;
	using DissociationReactionType = NEDissociationReaction;
	using ReSolutionReactionType = NEReSolutionReaction;
	using NucleationReactionType = NENucleationReaction;
	using SinkReactionType = NESinkReaction;

	using ReactionTypeList =
		std::tuple<ProductionReactionType, DissociationReactionType,
			ReSolutionReactionType, NucleationReactionType, SinkReactionType>;

	using ClusterGenerator = NEClusterGenerator;
	using ClusterUpdater = detail::NEClusterUpdater;
};

namespace detail
{
template <typename PlsmContext>
struct ClusterDataExtra<NEReactionNetwork, PlsmContext>
{
	using NetworkType = NEReactionNetwork;

	template <typename TData>
	using View = ViewType<TData, PlsmContext>;
	using MapType = Kokkos::UnorderedMap<int, int, PlsmContext>;

	using IndexType = detail::ReactionNetworkIndexType;

	ClusterDataExtra() = default;

	template <typename PC>
	KOKKOS_INLINE_FUNCTION
	ClusterDataExtra(const ClusterDataExtra<NetworkType, PC>& data) :
		constantRates(data.constantRates),
		fileClusterMap(data.fileClusterMap)
	{
	}

	template <typename PC>
	void
	deepCopy(const ClusterDataExtra<NetworkType, PC>& data)
	{
		if (!data.constantRates.is_allocated()) {
			return;
		}

		if (!constantRates.is_allocated()) {
			constantRates = create_mirror_view(data.constantRates);
		}

		if (!fileClusterMap.is_allocated()) {
			fileClusterMap = MapType(data.fileClusterMap.size());
		}

		deep_copy(constantRates, data.constantRates);
		deep_copy(fileClusterMap, data.fileClusterMap);
	}

	std::uint64_t
	getDeviceMemorySize() const noexcept
	{
		std::uint64_t ret = 0;

		ret += constantRates.required_allocation_size(constantRates.extent(0),
			constantRates.extent(1), constantRates.extent(2));
		ret += sizeof(fileClusterMap);

		return ret;
	}

	void
	initialize(IndexType size)
	{
		constantRates = View<double***>("Constant Rates", size, size + 1, 2);
		fileClusterMap = MapType(size);
	}

	View<double***> constantRates;
	MapType fileClusterMap;
};
} // namespace detail
} // namespace network
} // namespace core
} // namespace xolotl
