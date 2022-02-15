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
class TrapMutationReactionGenerator : public TBase
{
public:
	using Superclass = TBase;
	using NetworkType = typename TBase::NetworkType;
	using NetworkTraits = ReactionNetworkTraits<NetworkType>;
	using TrapMutationReactionType =
		typename NetworkTraits::TrapMutationReactionType;
	using IndexType = typename NetworkType::IndexType;
	using IndexView = typename Superclass::IndexView;
	using ClusterSetSubView = typename Superclass::ClusterSetSubView;
	using Count = typename Superclass::Count;
	using Construct = typename Superclass::Construct;

	TrapMutationReactionGenerator(const NetworkType& network);

	IndexType
	getRowMapAndTotalReactionCount();

	void
	setupCrsClusterSetSubView();

	KOKKOS_INLINE_FUNCTION
	void
	addTrapMutationReaction(Count, const ClusterSet& clusterSet) const;

	KOKKOS_INLINE_FUNCTION
	void
	addTrapMutationReaction(Construct, const ClusterSet& clusterSet) const;

	Kokkos::View<TrapMutationReactionType*>
	getTrapMutationReactions() const
	{
		return _tmReactions;
	}

	IndexType
	getNumberOfTrapMutationReactions() const
	{
		return _tmReactions.size();
	}

private:
	IndexView _clusterTMReactionCounts;

	IndexType _numPrecedingReactions{};
	IndexType _numTMReactions{};

	IndexView _tmCrsRowMap;
	ClusterSetSubView _tmCrsClusterSets;

	Kokkos::View<TrapMutationReactionType*> _tmReactions;
};

template <typename TNetwork, typename TReaction, typename TBase>
struct WrapTypeSpecificReactionGenerator<TNetwork, TReaction, TBase,
	std::enable_if_t<std::is_base_of_v<
		TrapMutationReaction<TNetwork, TReaction>, TReaction>>>
{
	using Type = TrapMutationReactionGenerator<TBase>;
};
} // namespace detail
} // namespace network
} // namespace core
} // namespace xolotl
