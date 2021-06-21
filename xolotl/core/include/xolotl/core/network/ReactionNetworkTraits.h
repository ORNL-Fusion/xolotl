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
/**
 * @brief Provide information about species and reaction types specific to a
 * network
 *
 * This is the place for things that must be specified explicitly for a given
 * network. Other template instantiations are specified in
 * detail::ReactionNetworkTypes and detail::ReactionNetworkProperties
 * based on the items specified here.
 * The default instantiation is empty. Therefore, it must be specialized for
 * every network type.
 *
 * The following items must be provided:
 * - **Species** (type)
 * - **numSpecies** (static constexpr std::size_t)
 * - **ProductionReactionType** (type)
 * - **DissociationReactionType** (type)
 * - (Optionally) **ReactionTypeList** (type; std::tuple of reaction types) if
 *   there are more than production and dissociation
 * - **ClusterGenerator** (type)
 * - (Optionally) **ClusterUpdater** (type)
 */
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

/**
 * @brief Customization mechanism for adding arbitrary network-specific data
 * items to ClusterData (ClusterDataImpl).
 *
 * The default is empty.
 *
 * The data here should be accessed only from code that "knows" what it
 * contains, that is, other network-specific code
 */
template <typename TNetwork, typename PlsmContext,
	template <typename> typename ViewConvert>
struct ClusterDataExtra;

template <typename TNetwork, typename PlsmContext>
struct ClusterDataHelper;

template <typename TNetwork, typename PlsmContext>
using ClusterData = typename ClusterDataHelper<TNetwork, PlsmContext>::Type;

template <typename TNetwork, typename PlsmContext>
struct ClusterDataRefHelper;

template <typename TNetwork, typename PlsmContext>
using ClusterDataRef =
	typename ClusterDataRefHelper<TNetwork, PlsmContext>::Type;

template <typename TNetwork>
struct ReactionDataRef;

template <typename TNetwork>
class DefaultClusterUpdater;

template <typename TNetwork>
class ReactionCollection;

template <typename TImpl, typename = std::void_t<>>
struct ClusterUpdaterHelper
{
	using Type = DefaultClusterUpdater<TImpl>;
};

/**
 * If ReactionNetworkTraits provides a ClusterUpdater member, then use that
 */
template <typename TImpl>
struct ClusterUpdaterHelper<TImpl,
	std::void_t<typename ReactionNetworkTraits<TImpl>::ClusterUpdater>>
{
	using Type = typename ReactionNetworkTraits<TImpl>::ClusterUpdater;
};

/**
 * @brief Instantiate types for things needed by the network based on the
 * ReactionNetworkTraits
 */
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

/**
 * @brief Label some properties of the network phase space based on the
 * ReactionNetworkTraits
 */
template <typename TImpl>
struct ReactionNetworkProperties
{
	using Traits = ReactionNetworkTraits<TImpl>;
	using Species = typename Traits::Species;
	static constexpr std::size_t numSpecies = Traits::numSpecies;
	using SpeciesSequence = SpeciesEnumSequence<Species, numSpecies>;
	static constexpr std::size_t numSpeciesNoI = SpeciesSequence::sizeNoI();
};

/**
 * By default, list production and dissociation types from ReactionNetworkTraits
 */
template <typename TNetwork, typename = std::void_t<>>
struct ReactionTypeListHelper
{
	using Traits = ReactionNetworkTraits<TNetwork>;
	using Type = std::tuple<typename Traits::ProductionReactionType,
		typename Traits::DissociationReactionType>;
};

/**
 * If ReactionNetworkTraits provides a ReactionTypeList member, then use that
 * See ReactionTypeListHelper
 */
template <typename TNetwork>
struct ReactionTypeListHelper<TNetwork,
	std::void_t<typename ReactionNetworkTraits<TNetwork>::ReactionTypeList>>
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
