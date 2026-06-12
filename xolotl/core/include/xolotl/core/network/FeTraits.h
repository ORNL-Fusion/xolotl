#pragma once

#include <xolotl/core/network/ReactionNetworkTraits.h>
#include <xolotl/core/network/detail/ClusterData.h>

namespace xolotl
{
namespace core
{
namespace network
{
class FeReactionNetwork;
class FeProductionReaction;
class FeDissociationReaction;
class FeSinkReaction;
class FeTrapReaction;
class FeClusterGenerator;

enum class FeSpeciesList
{
	Trap,
	He,
	V,
	I
};

inline const std::string&
toLabelString(FeSpeciesList species)
{
	static const std::string labelArray[] = {"Trap", "He", "V", "I"};
	return labelArray[static_cast<int>(species)];
}

inline const std::string&
toNameString(FeSpeciesList species)
{
	static const std::string nameArray[] = {
		"Trap", "Helium", "Vacancy", "Interstitial"};
	return nameArray[static_cast<int>(species)];
}

template <>
struct NumberOfSpecies<FeSpeciesList> : std::integral_constant<std::size_t, 4>
{
};

template <>
struct NumberOfInterstitialSpecies<FeSpeciesList> :
	std::integral_constant<std::size_t, 1>
{
};

template <>
struct NumberOfVacancySpecies<FeSpeciesList> :
	std::integral_constant<std::size_t, 1>
{
};

template <>
struct NumberOfTrapSpecies<FeSpeciesList> :
	std::integral_constant<std::size_t, 1>
{
};

template <>
struct SpeciesForGrouping<FeSpeciesList, 3>
{
	using Sequence = EnumSequence<FeSpeciesList, 3>;
	static constexpr auto first = Sequence(FeSpeciesList::He);
	static constexpr auto last = Sequence(FeSpeciesList::I);

	KOKKOS_INLINE_FUNCTION
	static constexpr std::underlying_type_t<FeSpeciesList>
	mapToMomentId(EnumSequence<FeSpeciesList, 3> value)
	{
		if (value == FeSpeciesList::I)
			return 2;
		return value();
	}
};

template <>
struct ReactionNetworkTraits<FeReactionNetwork>
{
	using Species = FeSpeciesList;

	static constexpr std::size_t numSpecies = 4;

	using ProductionReactionType = FeProductionReaction;
	using DissociationReactionType = FeDissociationReaction;
	using SinkReactionType = FeSinkReaction;
	using TrapReactionType = FeTrapReaction;

	using ReactionTypeList = std::tuple<ProductionReactionType,
		DissociationReactionType, SinkReactionType, TrapReactionType>;

	using ClusterGenerator = FeClusterGenerator;
};

namespace detail
{
template <typename PlsmContext>
struct ClusterDataExtra<FeReactionNetwork, PlsmContext>
{
	using NetworkType = FeReactionNetwork;

	template <typename TData>
	using View = ViewType<TData, PlsmContext>;

	using IndexType = detail::ReactionNetworkIndexType;

	ClusterDataExtra() = default;

	template <typename PC>
	KOKKOS_INLINE_FUNCTION
	ClusterDataExtra(const ClusterDataExtra<NetworkType, PC>& data) :
		dislocationCaptureRadius(data.dislocationCaptureRadius)
	{
	}

	template <typename PC>
	void
	deepCopy(const ClusterDataExtra<NetworkType, PC>& data)
	{
		if (!data.dislocationCaptureRadius.is_allocated()) {
			return;
		}

		if (!dislocationCaptureRadius.is_allocated()) {
			dislocationCaptureRadius =
				create_mirror_view(data.dislocationCaptureRadius);
		}

		deep_copy(dislocationCaptureRadius, data.dislocationCaptureRadius);
	}

	std::uint64_t
	getDeviceMemorySize() const noexcept
	{
		std::uint64_t ret = 0;

		ret += dislocationCaptureRadius.required_allocation_size(
			dislocationCaptureRadius.extent(0),
			dislocationCaptureRadius.extent(1));

		return ret;
	}

	void
	initialize(IndexType numClusters, IndexType gridSize = 0)
	{
		dislocationCaptureRadius =
			View<double**>("Dislocation Capture Radius", numClusters, 2);
	}

	View<double**> dislocationCaptureRadius;
};
} // namespace detail
} // namespace network
} // namespace core
} // namespace xolotl
