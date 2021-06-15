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
class BurstingReactionGenerator : public TBase
{
public:
	using Superclass = TBase;
	using NetworkType = typename TBase::NetworkType;
	using NetworkTraits = ReactionNetworkTraits<NetworkType>;
	using BurstingReactionType = typename NetworkTraits::BurstingReactionType;
	using IndexType = typename NetworkType::IndexType;
	using IndexView = typename Superclass::IndexView;
	using ClusterSetSubView = typename Superclass::ClusterSetSubView;
	using Count = typename Superclass::Count;
	using Construct = typename Superclass::Construct;

	BurstingReactionGenerator(const NetworkType& network);

	IndexType
	getRowMapAndTotalReactionCount();

	void
	setupCrsClusterSetSubView();

	KOKKOS_INLINE_FUNCTION
	void
	addBurstingReaction(Count, const ClusterSet& clusterSet) const;

	KOKKOS_INLINE_FUNCTION
	void
	addBurstingReaction(Construct, const ClusterSet& clusterSet) const;

	Kokkos::View<BurstingReactionType*>
	getBurstingReactions() const
	{
		return _burstingReactions;
	}

	IndexType
	getNumberOfBurstingReactions() const
	{
		return _burstingReactions.size();
	}

private:
	IndexView _clusterBurstingReactionCounts;

	IndexType _numPrecedingReactions{};
	IndexType _numBurstingReactions{};

	IndexView _burstingCrsRowMap;
	ClusterSetSubView _burstingCrsClusterSets;

	Kokkos::View<BurstingReactionType*> _burstingReactions;
};

template <typename TNetwork, typename TReaction, typename TBase>
struct WrapTypeSpecificReactionGenerator<TNetwork, TReaction, TBase,
	std::enable_if_t<
		std::is_base_of_v<BurstingReaction<TNetwork, TReaction>, TReaction>>>
{
	using Type = BurstingReactionGenerator<TBase>;
};
} // namespace detail
} // namespace network
} // namespace core
} // namespace xolotl
