#pragma once

#include <plsm/Subpaving.h>

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
template <typename TNetwork, typename PlsmContext>
class ClusterData;

template <typename TNetwork, typename PlsmContext>
class ClusterDataRef;

template <typename TImpl>
struct ReactionNetworkTypes
{
    using AmountType = typename std::uint32_t;
    using Traits = ReactionNetworkTraits<TImpl>;
    using Species = typename Traits::Species;
    using Subpaving = plsm::Subpaving<AmountType, Traits::numSpecies, Species>;
    using Region = typename Subpaving::RegionType;
    using Composition = typename Subpaving::PointType;
    using ClusterData = detail::ClusterData<TImpl, plsm::OnDevice>;
    using ClusterDataMirror = detail::ClusterData<TImpl, plsm::OnHost>;
    using ClusterDataRef = detail::ClusterDataRef<TImpl, plsm::OnDevice>;
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
}
}
}
