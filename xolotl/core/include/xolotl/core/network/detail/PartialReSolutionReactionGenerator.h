#pragma once

#include <xolotl/core/network/detail/ReactionGenerator.h>

namespace xolotl
{
namespace core
{
namespace network
{
namespace detail
{
/**
 * @brief Inherits from ReactionGeneratorBase and additionally fills collections
 * re-solution reactions depending on the subpaving.
 *
 * @tparam TBase The templated base class.
 */
template <typename TBase>
class PartialReSolutionReactionGenerator : public TBase
{
public:
	using Superclass = TBase;
	using NetworkType = typename TBase::NetworkType;
	using NetworkTraits = ReactionNetworkTraits<NetworkType>;
	using PartialReSolutionReactionType =
		typename NetworkTraits::PartialReSolutionReactionType;
	using IndexType = typename NetworkType::IndexType;
	using IndexView = typename Superclass::IndexView;
	using ClusterSetSubView = typename Superclass::ClusterSetSubView;
	using Count = typename Superclass::Count;
	using Construct = typename Superclass::Construct;

	PartialReSolutionReactionGenerator(const NetworkType& network);

	IndexType
	getRowMapAndTotalReactionCount();

	void
	setupCrsClusterSetSubView();

	KOKKOS_INLINE_FUNCTION
	void
	addPartialReSolutionReaction(Count, const ClusterSet& clusterSet) const;

	KOKKOS_INLINE_FUNCTION
	void
	addPartialReSolutionReaction(Construct, const ClusterSet& clusterSet) const;

	Kokkos::View<PartialReSolutionReactionType*>
	getPartialReSolutionReactions() const
	{
		return _reSoReactions;
	}

	IndexType
	getNumberOfPartialReSolutionReactions() const
	{
		return _reSoReactions.size();
	}

private:
	IndexView _clusterPartialReSoReactionCounts;

	IndexType _numPrecedingReactions{};
	IndexType _numPartialReSoReactions{};

	IndexView _reSoCrsRowMap;
	ClusterSetSubView _reSoCrsClusterSets;

	Kokkos::View<PartialReSolutionReactionType*> _reSoReactions;
};

template <typename TNetwork, typename TReaction, typename TBase>
struct WrapTypeSpecificReactionGenerator<TNetwork, TReaction, TBase,
	std::enable_if_t<std::is_base_of_v<
		PartialReSolutionReaction<TNetwork, TReaction>, TReaction>>>
{
	using Type = PartialReSolutionReactionGenerator<TBase>;
};
} // namespace detail
} // namespace network
} // namespace core
} // namespace xolotl
