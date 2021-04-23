#pragma once

#include <cstddef>
#include <cstdint>
#include <tuple>

#include <plsm/Subpaving.h>

#include <xolotl/config.h>
#include <xolotl/core/network/SpeciesEnumSequence.h>

namespace xolotl
{
namespace core
{
namespace network
{
template <typename TImpl>
struct ReactionNetworkTraits
{
};

namespace detail
{
using ReactionNetworkIndexType = ::xolotl::IdType;

using CompositionAmountType = std::uint32_t;

inline constexpr auto invalidNetworkIndex =
	plsm::invalid<ReactionNetworkIndexType>;

inline constexpr auto invalidSpeciesAmount =
	plsm::invalid<CompositionAmountType>;

template <typename TNetwork, typename PlsmContext>
struct ClusterDataHelper;

template <typename TNetwork, typename PlsmContext>
using ClusterData = typename ClusterDataHelper<TNetwork, PlsmContext>::Type;

template <typename TNetwork, typename PlsmContext>
struct ClusterDataRefHelper;

template <typename TNetwork, typename PlsmContext>
using ClusterDataRef =
	typename ClusterDataRefHelper<TNetwork, PlsmContext>::Type;

template <typename TNetwork, typename PlsmContext,
	template <typename> typename ViewConvert>
struct ClusterDataExtra
{
	ClusterDataExtra() = default;

	template <typename TOtherPlsmContext,
		template <typename> typename TOtherViewConvert>
	KOKKOS_INLINE_FUNCTION
	ClusterDataExtra(
		const ClusterDataExtra<TNetwork, TOtherPlsmContext, TOtherViewConvert>&)
	{
	}

	std::uint64_t
	getDeviceMemorySize() const noexcept
	{
		return 0;
	}
};

template <typename TN1, typename PC1, template <typename> typename VC1,
	typename TN2, typename PC2, template <typename> typename VC2>
inline void
deepCopy([[maybe_unused]] ClusterDataExtra<TN1, PC1, VC1> to,
	[[maybe_unused]] ClusterDataExtra<TN2, PC2, VC2> from)
{
}

template <typename TNetwork>
struct ReactionDataRef;

template <typename TNetwork>
class DefaultClusterUpdater;

template <typename TNetwork>
class ReactionCollection;

/*!
 * Stand-in for C++17 std::void_t
 */
template <typename...>
using VoidType = void;

template <typename TImpl, typename = VoidType<>>
struct ClusterUpdaterHelper
{
	using Type = DefaultClusterUpdater<TImpl>;
};

template <typename TImpl>
struct ClusterUpdaterHelper<TImpl,
	VoidType<typename ReactionNetworkTraits<TImpl>::ClusterUpdater>>
{
	using Type = typename ReactionNetworkTraits<TImpl>::ClusterUpdater;
};

template <typename TImpl>
struct ReactionNetworkTypes
{
	using IndexType = ReactionNetworkIndexType;
	using AmountType = CompositionAmountType;
	using Traits = ReactionNetworkTraits<TImpl>;
	using Species = typename Traits::Species;
	using Subpaving = plsm::Subpaving<AmountType, Traits::numSpecies, Species>;
	using Region = typename Subpaving::RegionType;
	using Composition = typename Subpaving::PointType;
	using ClusterData = detail::ClusterData<TImpl, plsm::OnDevice>;
	using ClusterDataMirror = detail::ClusterData<TImpl, plsm::OnHost>;
	using ClusterDataRef = detail::ClusterDataRef<TImpl, plsm::OnDevice>;
	using ClusterUpdater = typename ClusterUpdaterHelper<TImpl>::Type;
	using ReactionCollection = detail::ReactionCollection<TImpl>;
	using ReactionDataRef = typename detail::ReactionDataRef<TImpl>;
};

template <typename TImpl>
struct ReactionNetworkProperties
{
	using Traits = ReactionNetworkTraits<TImpl>;
	using Species = typename Traits::Species;
	static constexpr std::size_t numSpecies = Traits::numSpecies;
	using SpeciesSequence = SpeciesEnumSequence<Species, numSpecies>;
	static constexpr std::size_t numSpeciesNoI = SpeciesSequence::sizeNoI();
};

template <typename TNetwork, typename = VoidType<>>
struct ReactionTypeListHelper
{
	using Traits = ReactionNetworkTraits<TNetwork>;
	using Type = std::tuple<typename Traits::ProductionReactionType,
		typename Traits::DissociationReactionType>;
};

template <typename TNetwork>
struct ReactionTypeListHelper<TNetwork,
	VoidType<typename ReactionNetworkTraits<TNetwork>::ReactionTypeList>>
{
	using Type = typename ReactionNetworkTraits<TNetwork>::ReactionTypeList;
};

template <typename TNetwork>
using ReactionTypeList = typename ReactionTypeListHelper<TNetwork>::Type;

template <typename TSpeciesEnum>
KOKKOS_INLINE_FUNCTION
constexpr ReactionNetworkIndexType
toIndex(TSpeciesEnum species)
{
	static_assert(std::is_enum_v<TSpeciesEnum>);
	return static_cast<ReactionNetworkIndexType>(species);
}
} // namespace detail
} // namespace network
} // namespace core
} // namespace xolotl
