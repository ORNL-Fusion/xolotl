#pragma once

#include <tuple>

#include <xolotl/core/network/ReactionNetworkTraits.h>
#include <xolotl/core/network/detail/ClusterData.h>

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

namespace detail
{
template <typename PlsmContext>
struct ClusterDataExtra<AlloyReactionNetwork, PlsmContext>
{
	using NetworkType = AlloyReactionNetwork;

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
