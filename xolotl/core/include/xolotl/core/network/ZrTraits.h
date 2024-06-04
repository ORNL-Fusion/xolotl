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
class ZrReactionNetwork;
class ZrProductionReaction;
class ZrDissociationReaction;
class ZrSinkReaction;
class ZrConstantReaction;
class ZrClusterGenerator;
namespace detail
{
class ZrClusterUpdater;
}

enum class ZrSpecies
{
	V,
	Basal,
	I
};

inline const std::string&
toLabelString(ZrSpecies species)
{
	static const std::string labelArray[] = {"V", "Basal", "I"};
	return labelArray[static_cast<int>(species)];
}

inline const std::string&
toNameString(ZrSpecies species)
{
	static const std::string nameArray[] = {"Vacancy", "Basal", "Interstitial"};
	return nameArray[static_cast<int>(species)];
}

template <>
struct NumberOfInterstitialSpecies<ZrSpecies> :
	std::integral_constant<std::size_t, 1>
{
};

template <>
struct NumberOfVacancySpecies<ZrSpecies> :
	std::integral_constant<std::size_t, 2>
{
};

template <>
struct SpeciesForGrouping<ZrSpecies, 3>
{
	using Sequence = EnumSequence<ZrSpecies, 3>;
	static constexpr auto first = Sequence(ZrSpecies::V);
	static constexpr auto last = Sequence(ZrSpecies::I);

	KOKKOS_INLINE_FUNCTION
	static constexpr std::underlying_type_t<ZrSpecies>
	mapToMomentId(EnumSequence<ZrSpecies, 3>)
	{
		return 0;
	}
};

template <>
struct ReactionNetworkTraits<ZrReactionNetwork>
{
	using Species = ZrSpecies;

	static constexpr std::size_t numSpecies = 3;

	using ProductionReactionType = ZrProductionReaction;
	using DissociationReactionType = ZrDissociationReaction;
	using SinkReactionType = ZrSinkReaction;
	using ConstantReactionType = ZrConstantReaction;

	using ReactionTypeList = std::tuple<ProductionReactionType,
		DissociationReactionType, SinkReactionType, ConstantReactionType>;

	using ClusterGenerator = ZrClusterGenerator;
	using ClusterUpdater = detail::ZrClusterUpdater;
};

namespace detail
{
template <typename PlsmContext>
struct ClusterDataExtra<ZrReactionNetwork, PlsmContext>
{
	using NetworkType = ZrReactionNetwork;

	template <typename TData>
	using View = ViewType<TData, PlsmContext>;

	using IndexType = detail::ReactionNetworkIndexType;

	ClusterDataExtra() = default;

	template <typename PC>
	KOKKOS_INLINE_FUNCTION
	ClusterDataExtra(const ClusterDataExtra<NetworkType, PC>& data) :
		anisotropyRatio(data.anisotropyRatio),
		dislocationCaptureRadius(data.dislocationCaptureRadius)
	{
	}

	template <typename PC>
	void
	deepCopy(const ClusterDataExtra<NetworkType, PC>& data)
	{
		if (!data.anisotropyRatio.is_allocated()) {
			return;
		}

		if (!anisotropyRatio.is_allocated()) {
			anisotropyRatio = create_mirror_view(data.anisotropyRatio);
		}

		deep_copy(anisotropyRatio, data.anisotropyRatio);

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

		ret += anisotropyRatio.required_allocation_size(
			anisotropyRatio.extent(0), anisotropyRatio.extent(1));
		ret += dislocationCaptureRadius.required_allocation_size(
			dislocationCaptureRadius.extent(0),
			dislocationCaptureRadius.extent(1));

		return ret;
	}

	void
	initialize(IndexType numClusters, IndexType gridSize = 0)
	{
		anisotropyRatio =
			View<double**>("Anisotropy Ratio", numClusters, gridSize);
		dislocationCaptureRadius =
			View<double**>("Dislocation Capture Radius", numClusters, 2);
	}

	View<double**> anisotropyRatio;
	View<double**> dislocationCaptureRadius;
};
} // namespace detail
} // namespace network
} // namespace core
} // namespace xolotl
