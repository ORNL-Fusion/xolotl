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
template <typename TBase>
class ReSolutionReactionGenerator : public TBase
{
public:
	using Superclass = TBase;
	using NetworkType = typename TBase::NetworkType;
	using NetworkTraits = ReactionNetworkTraits<NetworkType>;
	using ReSolutionReactionType =
		typename NetworkTraits::ReSolutionReactionType;
	using IndexType = typename NetworkType::IndexType;
	using IndexView = typename Superclass::IndexView;
	using ClusterSetSubView = typename Superclass::ClusterSetSubView;
	using Count = typename Superclass::Count;
	using Construct = typename Superclass::Construct;

	ReSolutionReactionGenerator(const NetworkType& network);

	IndexType
	getRowMapAndTotalReactionCount();

	void
	setupCrsClusterSetSubView();

	KOKKOS_INLINE_FUNCTION
	void
	addReSolutionReaction(Count, const ClusterSet& clusterSet) const;

	KOKKOS_INLINE_FUNCTION
	void
	addReSolutionReaction(Construct, const ClusterSet& clusterSet) const;

	Kokkos::View<ReSolutionReactionType*>
	getReSolutionReactions() const
	{
		return _reSoReactions;
	}

	IndexType
	getNumberOfReSolutionReactions() const
	{
		return _reSoReactions.size();
	}

private:
	IndexView _clusterReSoReactionCounts;

	IndexType _numPrecedingReactions{};
	IndexType _numReSoReactions{};

	IndexView _reSoCrsRowMap;
	ClusterSetSubView _reSoCrsClusterSets;

	Kokkos::View<ReSolutionReactionType*> _reSoReactions;
};

template <typename TNetwork, typename TReaction, typename TBase>
struct WrapTypeSpecificReactionGenerator<TNetwork, TReaction, TBase,
	std::enable_if_t<std::is_base_of<ReSolutionReaction<TNetwork, TReaction>,
		TReaction>::value>>
{
	using Type = ReSolutionReactionGenerator<TBase>;
};
} // namespace detail
} // namespace network
} // namespace core
} // namespace xolotl
