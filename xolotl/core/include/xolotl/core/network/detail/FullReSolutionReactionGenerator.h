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
class FullReSolutionReactionGenerator : public TBase
{
public:
	using Superclass = TBase;
	using NetworkType = typename TBase::NetworkType;
	using NetworkTraits = ReactionNetworkTraits<NetworkType>;
	using FullReSolutionReactionType =
		typename NetworkTraits::FullReSolutionReactionType;
	using IndexType = typename NetworkType::IndexType;
	using IndexView = typename Superclass::IndexView;
	using ClusterSetSubView = typename Superclass::ClusterSetSubView;
	using Count = typename Superclass::Count;
	using Construct = typename Superclass::Construct;

	FullReSolutionReactionGenerator(const NetworkType& network);

	IndexType
	getRowMapAndTotalReactionCount();

	void
	setupCrsClusterSetSubView();

	KOKKOS_INLINE_FUNCTION
	void
	addFullReSolutionReaction(Count, const ClusterSet& clusterSet) const;

	KOKKOS_INLINE_FUNCTION
	void
	addFullReSolutionReaction(Construct, const ClusterSet& clusterSet) const;

	Kokkos::View<FullReSolutionReactionType*>
	getFullReSolutionReactions() const
	{
		return _reSoReactions;
	}

	IndexType
	getNumberOfFullReSolutionReactions() const
	{
		return _reSoReactions.size();
	}

private:
	IndexView _clusterFullReSoReactionCounts;

	IndexType _numPrecedingReactions{};
	IndexType _numFullReSoReactions{};

	IndexView _reSoCrsRowMap;
	ClusterSetSubView _reSoCrsClusterSets;

	Kokkos::View<FullReSolutionReactionType*> _reSoReactions;
};

template <typename TNetwork, typename TReaction, typename TBase>
struct WrapTypeSpecificReactionGenerator<TNetwork, TReaction, TBase,
	std::enable_if_t<std::is_base_of_v<
		FullReSolutionReaction<TNetwork, TReaction>, TReaction>>>
{
	using Type = FullReSolutionReactionGenerator<TBase>;
};
} // namespace detail
} // namespace network
} // namespace core
} // namespace xolotl
