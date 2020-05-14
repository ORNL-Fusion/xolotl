#pragma once

#include <cstddef>
#include <cstdint>
#include <tuple>

#include <plsm/Subpaving.h>

#include <experimental/ClusterData.h>
#include <experimental/SpeciesEnumSequence.h>

namespace xolotlCore
{
namespace experimental
{

template <typename TImpl>
struct ReactionNetworkTraits
{
};

namespace detail
{
using ReactionNetworkIndexType = std::size_t;

using CompositionAmountType = std::uint32_t;

template <typename TNetwork, typename PlsmContext>
class ClusterData;

template <typename TNetwork, typename PlsmContext>
class ClusterDataRef;

template <typename TImpl>
class DefaultClusterUpdater;

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

struct InvalidIndex
{
    static constexpr auto value = plsm::invalid<ReactionNetworkIndexType>;
};

struct InvalidSpeciesAmount
{
    static constexpr auto value = plsm::invalid<CompositionAmountType>;
};

template <typename TNetwork, typename = VoidType<>>
struct ReactionTypeListHelper
{
    using Traits = ReactionNetworkTraits<TNetwork>;
    using Type =
        std::tuple<typename Traits::ProductionReactionType,
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
}
}
}
